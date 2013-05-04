/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

/// Constructs this vector with the given values
SVECTOR6::SVECTOR6(REAL x, REAL y, REAL z, REAL a, REAL b, REAL c)
{
  _data[0] = x;
  _data[1] = y;
  _data[2] = z;
  _data[3] = a;
  _data[4] = b;
  _data[5] = c;
}

/// Constructs this vector with the given values
SVECTOR6::SVECTOR6(REAL x, REAL y, REAL z, REAL a, REAL b, REAL c, boost::shared_ptr<const POSE3> pose)
{
  _data[0] = x;
  _data[1] = y;
  _data[2] = z;
  _data[3] = a;
  _data[4] = b;
  _data[5] = c;
  this->pose = pose;
}

/// Constructs this vector from the given array
/**
 * \param array a 6-dimensional (or larger) array
 */
SVECTOR6::SVECTOR6(const REAL* array)
{
  for (unsigned i=0; i< 6; i++)
    _data[i] = array[i];
}

/// Constructs this vector from the given array
/**
 * \param array a 6-dimensional (or larger) array
 */
SVECTOR6::SVECTOR6(const REAL* array, boost::shared_ptr<const POSE3> pose)
{
  for (unsigned i=0; i< 6; i++)
    _data[i] = array[i];
  this->pose = pose;
}

/// Constructs the given spatial vector with given upper and lower components
SVECTOR6::SVECTOR6(const VECTOR3& upper, const VECTOR3& lower)
{
  set_upper(upper);
  set_lower(lower);
}

/// Constructs the given spatial vector with given upper and lower components
SVECTOR6::SVECTOR6(const VECTOR3& upper, const VECTOR3& lower, boost::shared_ptr<const POSE3> pose)
{
  set_upper(upper);
  set_lower(lower);
  this->pose = pose;
}

/// Gets an iterator to the beginning of the data
CONST_ITERATOR SVECTOR6::begin() const
{
  CONST_ITERATOR i;
  i._count = 0;
  i._sz = 6;
  i._ld = 6;
  i._rows = i._sz;
  i._columns = 1;
  i._data_start = i._current_data  = _data;
  return i;
}

/// Gets an iterator to the end of the data
CONST_ITERATOR SVECTOR6::end() const
{
  CONST_ITERATOR i;
  i._count = 6;
  i._sz = 6;
  i._ld = 6;
  i._rows = i._sz;
  i._columns = 1;
  i._data_start = _data;
  i._current_data  = data() + i._sz;
  return i;
}

/// Gets an iterator to the beginning of the data
ITERATOR SVECTOR6::begin()
{
  ITERATOR i;
  i._count = 0;
  i._sz = 6;
  i._ld = 6;
  i._rows = i._sz;
  i._columns = 1;
  i._data_start = i._current_data  = _data;
  return i;
}

/// Gets an iterator to the end of the data
ITERATOR SVECTOR6::end()
{
  ITERATOR i;
  i._count = 6;
  i._sz = 6;
  i._ld = 6;
  i._rows = i._sz;
  i._columns = 1;
  i._data_start = _data;
  i._current_data  = data() + i._sz;
  return i;
}

/// Computes the spatial cross product between two vectors
SVECTOR6 SVECTOR6::spatial_cross(const SVECTOR6& v1, const SVECTOR6& v2)
{
  // verify that both vectors are defined in the same frame
  #ifndef NEXCEPT
  if (v1.pose != v2.pose)
    throw FrameException();
  #endif

  VECTOR3 ax = v1.get_upper();
  VECTOR3 bx = v1.get_lower();

  // multiply
  VECTOR3 v2top = v2.get_upper();
  VECTOR3 v2bot = v2.get_lower();
  VECTOR3 top = VECTOR3::cross(ax, v2top);
  VECTOR3 bot = VECTOR3::cross(bx, v2top) + VECTOR3::cross(ax, v2bot);
  return SVECTOR6(top, bot);
}

/// Gets the lower 3-dimensional vector
VECTOR3 SVECTOR6::get_lower() const
{
  return VECTOR3(_data[3], _data[4], _data[5]);
}

/// Gets the upper 3-dimensional vector
VECTOR3 SVECTOR6::get_upper() const
{
  return VECTOR3(_data[0], _data[1], _data[2]);
}
 
/// Sets the lower 3-dimensional vector
void SVECTOR6::set_lower(const VECTOR3& lower)
{
  _data[3] = lower[0];
  _data[4] = lower[1];
  _data[5] = lower[2];
}

/// Sets the upper 3-dimensional vector
void SVECTOR6::set_upper(const VECTOR3& upper)
{
  _data[0] = upper[0];
  _data[1] = upper[1];
  _data[2] = upper[2];
}

/// Performs the dot product
REAL SVECTOR6::dot(const SVECTOR6& v1, const SVECTOR6& v2)
{
  #ifndef NEXCEPT
  if (v1.pose != v2.pose)
    throw FrameException();
  #endif

  return v1[3]*v2[0] + v1[4]*v2[1] + v1[5]*v2[2] + v1[0]*v2[3] + v1[1]*v2[4] + v1[2]*v2[5];
}

/// Copies this vector from another SVECTOR6
SVECTOR6& SVECTOR6::operator=(const SVECTOR6& v)
{
  pose = v.pose;
  _data[0] = v._data[0];
  _data[1] = v._data[1];
  _data[2] = v._data[2];
  _data[3] = v._data[3];
  _data[4] = v._data[4];
  _data[5] = v._data[5];
  return *this;
}

/// Returns the negation of this vector
SVECTOR6 SVECTOR6::operator-() const
{
  SVECTOR6 v;
  v._data[0] = -_data[0]; 
  v._data[1] = -_data[1]; 
  v._data[2] = -_data[2]; 
  v._data[3] = -_data[3]; 
  v._data[4] = -_data[4]; 
  v._data[5] = -_data[5]; 
  v.pose = pose;

  return v;
}

/// Multiplies this vector by a scalar in place
SVECTOR6& SVECTOR6::operator*=(REAL scalar)
{
  _data[0] *= scalar;
  _data[1] *= scalar;
  _data[2] *= scalar;
  _data[3] *= scalar;
  _data[4] *= scalar;
  _data[5] *= scalar;

  return *this;
}

/// Adds another vector to this one in place
SVECTOR6& SVECTOR6::operator+=(const SVECTOR6& v)
{
  #ifndef NEXCEPT
  if (pose != v.pose)
    throw FrameException();
  #endif

  _data[0] += v._data[0];
  _data[1] += v._data[1];
  _data[2] += v._data[2];
  _data[3] += v._data[3];
  _data[4] += v._data[4];
  _data[5] += v._data[5];

  return *this;
}

/// Subtracts another vector from this one in place
SVECTOR6& SVECTOR6::operator-=(const SVECTOR6& v)
{
  #ifndef NEXCEPT
  if (pose != v.pose)
    throw FrameException();
  #endif

  _data[0] -= v._data[0];
  _data[1] -= v._data[1];
  _data[2] -= v._data[2];
  _data[3] -= v._data[3];
  _data[4] -= v._data[4];
  _data[5] -= v._data[5];

  return *this;
}

/// Adds this vector to another and returns the result in a new vector
SVECTOR6 SVECTOR6::operator+(const SVECTOR6& v) const
{
  #ifndef NEXCEPT
  if (pose != v.pose)
    throw FrameException();
  #endif

  SVECTOR6 result;
  std::transform(begin(), end(), v.begin(), result.begin(), std::plus<REAL>());
  result.pose = pose;
  return result;
}

/// Subtracts another vector from this vector and returns the result in a new vector
SVECTOR6 SVECTOR6::operator-(const SVECTOR6& v) const
{
  #ifndef NEXCEPT
  if (pose != v.pose)
    throw FrameException();
  #endif

  SVECTOR6 result;
  std::transform(begin(), end(), v.begin(), result.begin(), std::minus<REAL>());
  result.pose = pose;
  return result;
}

