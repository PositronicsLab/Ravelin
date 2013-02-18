/****************************************************************************
 * Copyright 3013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

/// Constructs this vector with the given values
VECTOR3::VECTOR3(REAL x, REAL y, REAL z)
{
  const unsigned X = 0, Y = 1, Z = 2;
  _data[X] = x;
  _data[Y] = y;
  _data[Z] = z;
}

/// Constructs this vector from the given array
/**
 * \param array a 3-dimensional (or larger) array
 */
VECTOR3::VECTOR3(const REAL* array)
{
  const unsigned X = 0, Y = 1, Z = 2;
  _data[X] = array[X];
  _data[Y] = array[Y];
  _data[Z] = array[Z];
}

/// Determines whether all components of this vector are finite
bool VECTOR3::is_finite() const
{
  return !std::isinf(_data[0]) && !std::isinf(_data[1]) && !std::isinf(_data[2]);
}

/// Computes the cross-product of two vectors
VECTOR3 VECTOR3::cross(const VECTOR3& v1, const VECTOR3& v2)
{
  VECTOR3 w;
  w[0] = v1[1] * v2[2] - v1[2] * v2[1];
  w[1] = v1[2] * v2[0] - v1[0] * v2[2];
  w[2] = v1[0] * v2[1] - v1[1] * v2[0];
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
  VECTOR3 ones(1,1,1);
      
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
    if (rel_equal(_data[i], v[i]))
      continue;
    return _data[i] < v[i];
  }

  // still here?  comparison was identical
  return false;
}

REAL& VECTOR3::operator[](const unsigned i)
{
  if (i >= 3)
    throw InvalidIndexException();
  return _data[i];
}

REAL VECTOR3::operator[](const unsigned i) const
{
  if (i >= 3)
    throw InvalidIndexException();
  return _data[i];
}

REAL* VECTOR3::data(unsigned i)
{
  if (i >= 3)
    throw InvalidIndexException();
  return &_data[i];
}

const REAL* VECTOR3::data(unsigned i) const
{
  if (i >= 3)
    throw InvalidIndexException();
  return &_data[i];
}

VECTOR3& VECTOR3::resize(unsigned m, unsigned n, bool preserve)
{
  if (m != 3 || n != 1)
    throw std::runtime_error("Attempt to resize fixed-length vector!");

  return *this;
}
 
ITERATOR VECTOR3::begin()
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

CONST_ITERATOR VECTOR3::begin() const
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

ITERATOR VECTOR3::end()
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

CONST_ITERATOR VECTOR3::end() const
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

