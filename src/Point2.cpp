/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

/// Constructs this vector with the given values
POINT2::POINT2(REAL x, REAL y)
{
  const unsigned X = 0, Y = 1;
  _data[X] = x;
  _data[Y] = y;
}

/// Constructs this vector from the given origin 
POINT2::POINT2(const ORIGIN2& o)
{
  const unsigned X = 0, Y = 1;
  _data[X] = o.x();
  _data[Y] = o.y();
}

/// Constructs this vector with the given values
POINT2::POINT2(REAL x, REAL y, boost::shared_ptr<POSE2> p)
{
  const unsigned X = 0, Y = 1;
  _data[X] = x;
  _data[Y] = y;
  pose = p;
}

/// Constructs this vector from the given array
/**
 * \param array a 2-dimensional (or larger) array
 */
POINT2::POINT2(const REAL* array)
{
  const unsigned X = 0, Y = 1;
  _data[X] = array[X];
  _data[Y] = array[Y];
}

/// Constructs this point from the given origin
POINT2& POINT2::operator=(const ORIGIN2& o)
{
  x() = o.x();
  y() = o.y();
  return *this;
}

REAL& POINT2::operator[](const unsigned i)
{
  #ifndef NEXCEPT
  if (i >= 2)
    throw InvalidIndexException();
  #endif

  return _data[i];
}

const REAL& POINT2::operator[](const unsigned i) const
{
  #ifndef NEXCEPT
  if (i >= 2)
    throw InvalidIndexException();
  #endif

  return _data[i];
}

REAL* POINT2::data(unsigned i)
{
  #ifndef NEXCEPT
  if (i >= 2)
    throw InvalidIndexException();
  #endif

  return &_data[i];
}

const REAL* POINT2::data(unsigned i) const
{
  #ifndef NEXCEPT
  if (i >= 2)
    throw InvalidIndexException();
  #endif 

  return &_data[i];
}

POINT2& POINT2::operator=(const POINT2& p)
{
  x() = p.x();
  y() = p.y();
  pose = p.pose;
  return *this;
}

/// Adds a point and an origin to return a point
POINT2 POINT2::operator+(const ORIGIN2& o) const
{
  POINT2 result;
  result.x() = x() + o.x();
  result.y() = y() + o.y();
  result.pose = pose;
  return result;
}

/// Subtracts an origin from a point to return a point
POINT2 POINT2::operator-(const ORIGIN2& o) const
{
  POINT2 result;
  result.x() = x() - o.x();
  result.y() = y() - o.y();
  result.pose = pose;
  return result;
}

POINT2 POINT2::operator+(const POINT2& p) const
{
  #ifndef NEXCEPT
  if (pose != p.pose)
    throw FrameException();
  #endif
  POINT2 result;
  result.x() = x() + p.x();
  result.y() = y() + p.y();
  result.pose = p.pose;
  return result;
}

POINT2 POINT2::operator-(const POINT2& p) const
{
  #ifndef NEXCEPT
  if (pose != p.pose)
    throw FrameException();
  #endif
  POINT2 result;
  result.x() = x() - p.x();
  result.y() = y() - p.y();
  result.pose = p.pose;
  return result;
}

POINT2& POINT2::operator+=(const POINT2& p) 
{
  #ifndef NEXCEPT
  if (pose != p.pose)
    throw FrameException();
  #endif
  x() += p.x();
  y() += p.y();
  return *this;
}

/// Adds an origin to this
POINT2& POINT2::operator+=(const ORIGIN2& o) 
{
  x() += o.x();
  y() += o.y();
  return *this;
}

POINT2& POINT2::operator-=(const POINT2& p) 
{
  #ifndef NEXCEPT
  if (pose != p.pose)
    throw FrameException();
  #endif
  x() -= p.x();
  y() -= p.y();
  return *this;
}

/// Subtracts an origin from this
POINT2& POINT2::operator-=(const ORIGIN2& o) 
{
  x() -= o.x();
  y() -= o.y();
  return *this;
}

