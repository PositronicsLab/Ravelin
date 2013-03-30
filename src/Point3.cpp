/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

/// Constructs this vector with the given values
POINT3::POINT3(REAL x, REAL y, REAL z)
{
  const unsigned X = 0, Y = 1, Z = 2;
  _data[X] = x;
  _data[Y] = y;
  _data[Z] = z;
}

/// Constructs this vector with the given values and defined with respect to the given pose
POINT3::POINT3(REAL x, REAL y, REAL z, boost::shared_ptr<POSE3> p)
{
  const unsigned X = 0, Y = 1, Z = 2;
  _data[X] = x;
  _data[Y] = y;
  _data[Z] = z;
  pose = p;
}

/// Constructs this vector from the given array
/**
 * \param array a 3-dimensional (or larger) array
 */
POINT3::POINT3(const REAL* array)
{
  const unsigned X = 0, Y = 1, Z = 2;
  _data[X] = array[X];
  _data[Y] = array[Y];
  _data[Z] = array[Z];
}

/// Compares the two vectors lexographically
bool POINT3::operator<(const POINT3& v) const
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

REAL& POINT3::operator[](const unsigned i)
{
  #ifndef NEXCEPT
  if (i >= 3)
    throw InvalidIndexException();
  #endif
  return _data[i];
}

const REAL& POINT3::operator[](const unsigned i) const
{
  #ifndef NEXCEPT
  if (i >= 3)
    throw InvalidIndexException();
  #endif
  return _data[i];
}

REAL* POINT3::data(unsigned i)
{
  #ifndef NEXCEPT
  if (i >= 3)
    throw InvalidIndexException();
  #endif
  return &_data[i];
}

const REAL* POINT3::data(unsigned i) const
{
  #ifndef NEXCEPT
  if (i >= 3)
    throw InvalidIndexException();
  #endif
  return &_data[i];
}

POINT3& POINT3::operator=(const ORIGIN3& o)
{
  x() = o.x();
  y() = o.y();
  z() = o.z();
  return *this;
}

POINT3& POINT3::operator=(const POINT3& p)
{
  x() = p.x();
  y() = p.y();
  z() = p.z();
  pose = p.pose;
  return *this;
}

POINT3 POINT3::operator+(const POINT3& p) const
{
  #ifndef NEXCEPT
  if (pose != p.pose)
    throw FrameException();
  #endif
  POINT3 result;
  result.x() = x() + p.x();
  result.y() = y() + p.y();
  result.z() = z() + p.z();
  result.pose = p.pose;
  return result;
}

POINT3 POINT3::operator-(const POINT3& p) const
{
  #ifndef NEXCEPT
  if (pose != p.pose)
    throw FrameException();
  #endif
  POINT3 result;
  result.x() = x() - p.x();
  result.y() = y() - p.y();
  result.z() = z() - p.z();
  result.pose = p.pose;
  return result;
}

POINT3& POINT3::operator+=(const POINT3& p) 
{
  #ifndef NEXCEPT
  if (pose != p.pose)
    throw FrameException();
  #endif
  x() += p.x();
  y() += p.y();
  z() += p.z();
  return *this;
}

POINT3& POINT3::operator-=(const POINT3& p) 
{
  #ifndef NEXCEPT
  if (pose != p.pose)
    throw FrameException();
  #endif
  x() -= p.x();
  y() -= p.y();
  z() -= p.z();
  return *this;
}

