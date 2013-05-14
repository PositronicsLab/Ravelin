/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

/// Default constructor
/**
 * Sets matrix to the identity matrix
 */
TRANSFORM2::TRANSFORM2()
{
  set_identity();
}

/// Constructs a 2D pose from a rotation and translation vector
TRANSFORM2::TRANSFORM2(const ROT2& r, const ORIGIN2& v)
{
  set(r, v);
}

/// Constructs a 2D pose from a rotation and zero translation
TRANSFORM2::TRANSFORM2(const ROT2& r)
{
  set(r, ORIGIN2::zero());
}

/// Constructs a 2D pose using identity orientation and a translation vector
TRANSFORM2::TRANSFORM2(const ORIGIN2& v)
{
  set(ROT2::identity(), v);
}

TRANSFORM2& TRANSFORM2::operator=(const TRANSFORM2& p)
{
  r = p.r;
  x = p.x;
  source = p.source;
  target = p.target;
  return *this;
}

REAL TRANSFORM2::wrap(REAL theta)
{
  if (theta > (REAL) M_PI)
  {
    do
    {
      theta -= (REAL) M_PI;
    }
    while (theta > (REAL) M_PI);
  }
  else if (theta < (REAL) -M_PI)
  {
    do
    {
      theta += (REAL) M_PI;
    }
    while (theta < (REAL) -M_PI);
  }

  return theta;
}

/// Determines whether two poses in 2D are relatively equivalent
bool TRANSFORM2::rel_equal(const TRANSFORM2& p1, const TRANSFORM2& p2, REAL tol)
{
  // verify relative frames are identity 
  if (!(p1.source == p2.source && p1.target == p2.target))
    throw FrameException();

  // wrap two thetas to [-pi, pi]
  REAL theta1 = wrap(p1.r.theta);
  REAL theta2 = wrap(p2.r.theta);

  // check all components
  return OPS::rel_equal(p1.x[0], p2.x[0], tol) && 
         OPS::rel_equal(p1.x[1], p2.x[1], tol) &&
         OPS::rel_equal(theta1, theta2, tol);
}

/// Sets the pose from a rotation and translation vector
TRANSFORM2& TRANSFORM2::set(const ROT2& r, const ORIGIN2& v)
{
  this->r = r;
  this->x = v;
  return *this;
}

/// Sets this pose from a rotation only (translation will be zero'd) 
TRANSFORM2& TRANSFORM2::set(const ROT2& r)
{
  this->r = r;
  this->x = ORIGIN2::zero();
  return *this;
}

/// Sets this matrix to identity
TRANSFORM2& TRANSFORM2::set_identity()
{
  this->r = ROT2::identity();
  this->x = ORIGIN2::zero();
  return *this;
}

/// Transforms a vector from one pose to another 
VECTOR2 TRANSFORM2::transform(const VECTOR2& v) const
{
  #ifndef NEXCEPT
  if (v.pose != source)
    throw FrameException();
  #endif

  VECTOR2 result = r * v;
  result.pose = target;
  return result;
}

/// Transforms a vector from one pose to another 
VECTOR2 TRANSFORM2::inverse_transform(const VECTOR2& v) const
{
  #ifndef NEXCEPT
  if (v.pose != target)
    throw FrameException();
  #endif

  VECTOR2 result = ROT2::invert(r) * v;
  result.pose = source; 

  return result;
}

/// Transforms a point from one pose to another 
POINT2 TRANSFORM2::transform(const POINT2& p) const
{
  #ifndef NEXCEPT
  if (p.pose != source)
    throw FrameException();
  #endif

  POINT2 result = r * p + x;
  result.pose = target;
  return result;
}

/// Transforms a point from one pose to another 
POINT2 TRANSFORM2::inverse_transform(const POINT2& p) const
{
  #ifndef NEXCEPT
  if (p.pose != target)
    throw FrameException();
  #endif

  POINT2 result = ROT2::invert(r) * (p - x);
  result.pose = source;
  return result;
}

/// Special method for inverting a 2D pose 
TRANSFORM2 TRANSFORM2::inverse(const TRANSFORM2& p)
{
  return TRANSFORM2(p).invert();
}

/// Concatenates transformations
/**
 * If this transforms from frame b to frame c and T transforms from frame a
 * to frame b, then concatenation transforms from frame a to frame c.
 */ 
TRANSFORM2 TRANSFORM2::operator*(const TRANSFORM2& T) const
{
  #ifndef NEXCEPT
  if (source != T.target)
    throw FrameException();
  #endif
  TRANSFORM2 result(r * T.r, r*T.x + x);
  result.source = T.source;
  result.target = target;
  return result;
}

/// Outputs this matrix to the stream
std::ostream& Ravelin::operator<<(std::ostream& out, const TRANSFORM2& p)
{
  out << "r: " << p.r << " " << p.x << std::endl;
   
  return out;
}

