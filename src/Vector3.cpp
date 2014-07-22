/****************************************************************************
 * Copyright 3013 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

/// Constructs this vector with the given values
VECTOR3::VECTOR3(REAL x, REAL y, REAL z, boost::shared_ptr<const POSE3> pose)
{
  const unsigned X = 0, Y = 1, Z = 2;
  _data[X] = x;
  _data[Y] = y;
  _data[Z] = z;
  this->pose = pose;
}

/// Constructs this vector with the given values
VECTOR3::VECTOR3(REAL x, REAL y, REAL z, boost::shared_ptr<POSE3> pose)
{
  const unsigned X = 0, Y = 1, Z = 2;
  _data[X] = x;
  _data[Y] = y;
  _data[Z] = z;
  this->pose = boost::const_pointer_cast<const POSE3>(pose);
}

/// Constructs this vector from the given array
/**
 * \param array a 3-dimensional (or larger) array
 */
VECTOR3::VECTOR3(const REAL* array, boost::shared_ptr<const POSE3> pose)
{
  const unsigned X = 0, Y = 1, Z = 2;
  _data[X] = array[X];
  _data[Y] = array[Y];
  _data[Z] = array[Z];
  this->pose = pose;
}

/// Constructs this vector from the given array
/**
 * \param array a 3-dimensional (or larger) array
 */
VECTOR3::VECTOR3(const REAL* array, boost::shared_ptr<POSE3> pose)
{
  const unsigned X = 0, Y = 1, Z = 2;
  _data[X] = array[X];
  _data[Y] = array[Y];
  _data[Z] = array[Z];
  this->pose = boost::const_pointer_cast<const POSE3>(pose);
}

/// Determines whether all components of this vector are finite
bool VECTOR3::is_finite() const
{
  return !std::isinf(_data[0]) && !std::isinf(_data[1]) && !std::isinf(_data[2]);
}

/// Computes the cross-product of two vectors
VECTOR3 VECTOR3::cross(const VECTOR3& v1, const VECTOR3& v2)
{
  #ifndef NEXCEPT
  if (v1.pose != v2.pose)
    throw FrameException();
  #endif
  VECTOR3 w;
  w[0] = v1[1] * v2[2] - v1[2] * v2[1];
  w[1] = v1[2] * v2[0] - v1[0] * v2[2];
  w[2] = v1[0] * v2[1] - v1[1] * v2[0];
  w.pose = v1.pose;
  return w;
}

/// Determines a vector orthogonal to v
/**
 * \note the returned vector is not normalized
 */ 
VECTOR3 VECTOR3::determine_orthogonal_vec(const VECTOR3& v)
{
  const unsigned X = 0, Y = 1, Z = 2;
  
  // make a vector of all ones
  VECTOR3 ones(1,1,1,v.pose);
      
  // get the absolute values of the three components
  REAL x = std::fabs(v[X]);
  REAL y = std::fabs(v[Y]);
  REAL z = std::fabs(v[Z]);
    
  // make the component zero that is equal to the largest component of v1 
  if (x > y)
  {
    if (x > z)
      ones[X] = 0;
    else
      ones[Z] = 0;
  }
  else
  {
    if (y > z)
      ones[Y] = 0;
    else
      ones[Z] = 0;
  }
      
  // compute the cross product of v and the ones vector
  return VECTOR3::cross(v, ones);
}

/// Computes an orthonormal basis, given a single vector
/**
 * \return an orthonormal basis constructed from a single 3D vector, v1; v1 x v2 = v3
 */
void VECTOR3::determine_orthonormal_basis(const VECTOR3& v1, VECTOR3& v2, VECTOR3& v3)
{
  // get a vector orthogonal to v1
  v2 = VECTOR3::normalize(determine_orthogonal_vec(v1));
  
  // compute the second vector
  v3 = VECTOR3::normalize(VECTOR3::cross(v1, v2));
}

/// Compares the two vectors lexographically
bool VECTOR3::operator<(const VECTOR3& v) const
{
  const unsigned LEN = 3;
  for (unsigned i=0; i< LEN; i++)
  {
    if (OPS::rel_equal(_data[i], v[i]))
      continue;
    return _data[i] < v[i];
  }

  // still here?  comparison was identical
  return false;
}

REAL& VECTOR3::operator[](const unsigned i)
{
  #ifndef NEXCEPT
  if (i >= 3)
    throw InvalidIndexException();
  #endif
  return _data[i];
}

const REAL& VECTOR3::operator[](const unsigned i) const
{
  #ifndef NEXCEPT
  if (i >= 3)
    throw InvalidIndexException();
  #endif
  return _data[i];
}

REAL* VECTOR3::data(unsigned i)
{
  #ifndef NEXCEPT
  if (i >= 3)
    throw InvalidIndexException();
  #endif
  return &_data[i];
}

const REAL* VECTOR3::data(unsigned i) const
{
  #ifndef NEXCEPT
  if (i >= 3)
    throw InvalidIndexException();
  #endif
  return &_data[i];
}

VECTOR3& VECTOR3::resize(unsigned m, unsigned n, bool preserve)
{
  #ifndef NEXCEPT
  if (!((m == 3 && n == 1) || (m == 1 && n == 3)))
    throw std::runtime_error("Attempt to resize fixed-length vector!");
  #endif

  return *this;
}

/// Computes the dot product between two vectors
REAL VECTOR3::dot(const VECTOR3& v1, const VECTOR3& v2)
{
  #ifndef NEXCEPT
  if (v1.pose != v2.pose)
    throw FrameException();
  #endif
  return v1.x()*v2.x() + v1.y()*v2.y() + v1.z()*v2.z();
}

/// Adds two vectors in the same frame
VECTOR3& VECTOR3::operator+=(const VECTOR3& v)
{
  #ifndef NEXCEPT
  if (pose != v.pose)
    throw FrameException();
  #endif

  x() += v.x();
  y() += v.y();
  z() += v.z();
  return *this;
}
 
/// Adds a vector and an origin 
VECTOR3& VECTOR3::operator+=(const ORIGIN3& o)
{
  x() += o.x();
  y() += o.y();
  z() += o.z();
  return *this;
}

 /// Subtracts a vector from this (in the same frame)
VECTOR3& VECTOR3::operator-=(const VECTOR3& v)
{
  #ifndef NEXCEPT
  if (pose != v.pose)
    throw FrameException();
  #endif

  x() -= v.x();
  y() -= v.y();
  z() -= v.z();
  return *this;
}
 
/// Subtracts an origin from this vector
VECTOR3& VECTOR3::operator-=(const ORIGIN3& o)
{
  x() -= o.x();
  y() -= o.y();
  z() -= o.z();
  return *this;
}

COLUMN_ITERATOR VECTOR3::column_iterator_begin()
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

CONST_COLUMN_ITERATOR VECTOR3::column_iterator_begin() const
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

COLUMN_ITERATOR VECTOR3::column_iterator_end()
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

CONST_COLUMN_ITERATOR VECTOR3::column_iterator_end() const
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

ROW_ITERATOR VECTOR3::row_iterator_begin()
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

CONST_ROW_ITERATOR VECTOR3::row_iterator_begin() const
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

ROW_ITERATOR VECTOR3::row_iterator_end()
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

CONST_ROW_ITERATOR VECTOR3::row_iterator_end() const
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

