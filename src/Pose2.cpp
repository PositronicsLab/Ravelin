/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

/// Default constructor
/**
 * Sets matrix to the identity matrix
 */
POSE2::POSE2()
{
  set_identity();
}

/// Constructs a 2D pose from a rotation and translation vector
POSE2::POSE2(const ROT2& r, const ORIGIN2& v)
{
  set(r, v);
}

/// Constructs a 2D pose from a rotation and zero translation
POSE2::POSE2(const ROT2& r)
{
  set(r, ORIGIN2::zero());
}

/// Constructs a 2D pose using identity orientation and a translation vector
POSE2::POSE2(const ORIGIN2& v)
{
  set(ROT2::identity(), v);
}

POSE2& POSE2::operator=(const POSE2& p)
{
  r = p.r;
  x = p.x;
  return *this;
}

/// Interpolates between two 4x4 transforms using spherical linear interpolation
/**
 * \param T1 the matrix to use when t=0
 * \param T2 the matrix to use when t=1
 * \param t a real value in the interval [0,1]
 * \return the interpolated transform
 */
POSE2 POSE2::interpolate(const POSE2& T1, const POSE2& T2, REAL t)
{
  // interpolate the positions
  ORIGIN2 x = T1.x*(1-t) + T2.x*t;

  // interpolate the rotations
  ROT2 r = T1.r.theta*(1-t) + T2.r.theta*t;

  // return the new matrix
  return POSE2(r, x);
}

/// Sets the pose from a rotation and translation vector
POSE2& POSE2::set(const ROT2& r, const ORIGIN2& v)
{
  this->r = r;
  this->x = v;
  return *this;
}

/// Sets this pose from a rotation only (translation will be zero'd) 
POSE2& POSE2::set(const ROT2& r)
{
  this->r = r;
  this->x = ORIGIN2::zero();
  return *this;
}

/// Sets this matrix to identity
POSE2& POSE2::set_identity()
{
  this->r = ROT2::identity();
  this->x = ORIGIN2::zero();
  return *this;
}

/// Transforms a vector from one pose to another 
VECTOR2 POSE2::transform(boost::shared_ptr<const POSE2> p, const VECTOR2& v) const
{
  return transform(shared_from_this(), p, v);
}

/// Applies this pose to a vector 
VECTOR2 POSE2::transform(boost::shared_ptr<const POSE2> source, boost::shared_ptr<const POSE2> target, const VECTOR2& v) 
{
  #ifndef NEXCEPT
  if (source != v.pose)
    throw FrameException();
  #endif

  // compute the relative transform
  std::pair<ROT2, ORIGIN2> Tx = calc_transform(source, target);

  return Tx.first * v;
}

/// Transforms a point from one pose to another 
POINT2 POSE2::transform(boost::shared_ptr<const POSE2> p, const POINT2& point) const
{
  return transform(shared_from_this(), p, point);
}

/// Transforms a point from one pose to another 
POINT2 POSE2::transform(boost::shared_ptr<const POSE2> source, boost::shared_ptr<const POSE2> target, const POINT2& point)
{
  #ifndef NEXCEPT
  if (source != point.pose)
    throw FrameException();
  #endif

  // compute the relative transform
  std::pair<ROT2, ORIGIN2> Tx = calc_transform(source, target);

  // do the transform
  return Tx.first * point + Tx.second;
}

/// Special method for inverting a 2D pose in place
POSE2& POSE2::invert()
{
  // get the new rotation 
  r.inverse();

  // determine the new translation
  x = r * (-x);

  return *this;
}

/// Special method for inverseing a 4x4 transformation matrix
POSE2 POSE2::inverse(const POSE2& p)
{
  return POSE2(p).invert();
}

/// Transforms this pose by another
POSE2 POSE2::operator*(const POSE2& p) const
{
  #ifndef NEXCEPT
  if (rpose != p.rpose)
    throw FrameException();
  #endif
  return POSE2(r * p.r, r*p.x + x);
}

/// Determines whether pose p exists in the chain of relative poses (and, if so, how far down the chain)
bool POSE2::is_common(boost::shared_ptr<const POSE2> x, boost::shared_ptr<const POSE2> p, unsigned& i)
{
  // reset i
  i = 0;

  while (x)
  {
    if (x == p)
      return true;
    x = x->rpose;
    i++;
  }

  return false;
}

/// Computes the relative transformation from this pose to another
std::pair<ROT2, ORIGIN2> POSE2::calc_transform(boost::shared_ptr<const POSE2> p) const
{
  return calc_transform(shared_from_this(), p);
}

/// Computes the relative transformation from this pose to another
std::pair<ROT2, ORIGIN2> POSE2::calc_transform(boost::shared_ptr<const POSE2> source, boost::shared_ptr<const POSE2> target)
{
  std::pair<ROT2, ORIGIN2> result;
  boost::shared_ptr<const POSE2> r, s; 

  // check for special case: transformation to and from global frame
  if (source == target && !source)
  {
    result.first = ROT2::identity();
    result.second.set_zero();
    return result;
  }

  // check for special case: transformation to global frame
  if (!target)
  {
    // combine transforms from this to i: this will give aTl
    ROT2 left_r = source->r;
    ORIGIN2 left_x = source->x;
    s = source;
    while (s)
    {
      s = s->rpose;
      if (!s)
        break;
      left_x = s->x + s->r * left_x;
      left_r = s->r * left_r;
    }

    // setup the transform 
    result.first = left_r;      
    result.second = -left_x;
  }

  // check for special case: transformation from global frame
  if (!source)
  {
    // combine transforms from target to r
    ROT2 right_r = target->r;
    ORIGIN2 right_x = target->x;
    r = target;
    while (r)
    {
      r = r->rpose;
      if (!r)
        break;
      right_x = r->x + r->r * right_x; 
      right_r = r->r * right_r;
    }

    // compute the inverse pose of the right 
    ROT2 inv_right_r = ROT2::invert(right_r);
    ORIGIN2 inv_right_x = inv_right_r * (-right_x);

    // multiply the inverse pose of p by this 
    result.first = inv_right_r;      
    result.second = inv_right_r * inv_right_x;
    return result;
  }

  // if both transforms are defined relative to the same frame, this is easy
  if (source->rpose == target->rpose)
  {
    // compute the inverse pose of p 
    ROT2 target_r = ROT2::invert(target->r);
    ORIGIN2 target_x = target_r * (-target->x);

    // multiply the inverse pose of p by this 
    result.first = target_r * source->r;      
    result.second = target_r * (target_x - source->x);
  }
  else
  {
    // search for the common link
    unsigned i = std::numeric_limits<unsigned>::max();
    r = target;
    while (true)
    {
      if (is_common(source, r, i))
        break;
      else
      {
        assert(r);
        r = r->rpose;
      } 
    } 
    
     // combine transforms from this to i: this will give aTl
    ROT2 left_r = source->r;
    ORIGIN2 left_x = source->x;
    s = source;
    for (unsigned j=0; j < i; j++)
    {
      s = s->rpose;
      left_x = s->x + s->r * left_x;
      left_r = s->r * left_r;
    }

    // combine transforms from target to r
    ROT2 right_r = target->r;
    ORIGIN2 right_x = target->x;
    while (target != r)
    {
      target = target->rpose;
      right_x = target->x + target->r * right_x; 
      right_r = target->r * right_r;
    }

    // compute the inverse pose of the right 
    ROT2 inv_right_r = ROT2::invert(right_r);
    ORIGIN2 inv_right_x = inv_right_r * (-right_x);

    // multiply the inverse pose of p by this 
    result.first = inv_right_r * left_r;      
    result.second = inv_right_r * (inv_right_x - left_x);
  }

  return result;
}

/// Outputs this matrix to the stream
std::ostream& Ravelin::operator<<(std::ostream& out, const POSE2& p)
{
  out << "r: " << p.r << " " << p.x << std::endl;
   
  return out;
}

