/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

/// Default constructor
/**
 * Sets matrix to the identity matrix
 */
POSE2::POSE2(boost::shared_ptr<const POSE2> relative_pose)
{
  set_identity();
  rpose = relative_pose;
}

/// Constructs a 2D pose from a rotation and translation vector
POSE2::POSE2(const ROT2& r, const ORIGIN2& v, boost::shared_ptr<const POSE2> relative_pose)
{
  set(r, v);
  rpose = relative_pose;
}

/// Constructs a 2D pose from a rotation and zero translation
POSE2::POSE2(const ROT2& r, boost::shared_ptr<const POSE2> relative_pose)
{
  set(r, ORIGIN2::zero());
  rpose = relative_pose;
}

/// Constructs a 2D pose using identity orientation and a translation vector
POSE2::POSE2(const ORIGIN2& v, boost::shared_ptr<const POSE2> relative_pose)
{
  set(ROT2::identity(), v);
  rpose = relative_pose;
}

POSE2& POSE2::operator=(const POSE2& p)
{
  r = p.r;
  x = p.x;
  rpose = p.rpose;
  return *this;
}

/// Wraps an angle to [-pi, pi]
REAL POSE2::wrap(REAL theta)
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
bool POSE2::rel_equal(const POSE2& p1, const POSE2& p2, REAL tol)
{
  // verify relative frames are identity 
  if (p1.rpose != p2.rpose)
    throw FrameException();

  // wrap two thetas to [-pi, pi]
  REAL theta1 = wrap(p1.r.theta);
  REAL theta2 = wrap(p2.r.theta);

  // check all components
  return OPS::rel_equal(p1.x[0], p2.x[0], tol) && 
         OPS::rel_equal(p1.x[1], p2.x[1], tol) &&
         OPS::rel_equal(theta1, theta2, tol);
}

/// Tranforms a vector with an interpolated pose (between poses P1 and P2)
/**
 * \param P1 the pose to use when t=0
 * \param P2 the pose to use when t=1
 * \param t interpolation value
 * \param o the vector to transform 
 * \return the transformed vector 
 */
VECTOR2 POSE2::interpolate_transform_vector(const POSE2& P1, const POSE2& P2, REAL t, const ORIGIN2& o)
{
  #ifndef NEXCEPT
  if (P1.rpose != P2.rpose)
    throw FrameException();
  #endif

  // wrap two thetas to [-pi, pi]
  REAL theta1 = wrap(P1.r.theta);
  REAL theta2 = wrap(P2.r.theta);

  // interpolate the rotations
  ROT2 r = P1.r.theta*(1-t) + P2.r.theta*t;

  // transform the vector
  return VECTOR2(r*o, P1.rpose);
}

/// Tranforms a point with an interpolated pose (between poses P1 and P2)
/**
 * \param P1 the pose to use when t=0
 * \param P2 the pose to use when t=1
 * \param t interpolation value
 * \param o the point to transform 
 * \return the transformed point 
 */
VECTOR2 POSE2::interpolate_transform_point(const POSE2& P1, const POSE2& P2, REAL t, const ORIGIN2& o)
{
  #ifndef NEXCEPT
  if (P1.rpose != P2.rpose)
    throw FrameException();
  #endif

  // interpolate the positions
  ORIGIN2 x = P1.x*(1-t) + P2.x*t;

  // wrap two thetas to [-pi, pi]
  REAL theta1 = wrap(P1.r.theta);
  REAL theta2 = wrap(P2.r.theta);

  // interpolate the rotations
  ROT2 r = P1.r.theta*(1-t) + P2.r.theta*t;

  // trnsform the point 
  return VECTOR2(x+r*o, P1.rpose);
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
VECTOR2 POSE2::transform_vector(const VECTOR2& v) const
{
  #ifndef NEXCEPT
  boost::shared_ptr<const POSE2> pose;
  try
  {
    pose = shared_from_this();
  }
  catch (boost::bad_weak_ptr e)
  {
    std::cerr << "Pose2 can't check frames when allocated on the stack! Allocate on heap or disable debugging exceptions." << std::endl;
  }
  if (v.pose != pose)
    throw FrameException();
  #endif

  VECTOR2 result = r * v;
  result.pose = rpose;
  return result;
}

/// Transforms a vector from one pose to another 
VECTOR2 POSE2::inverse_transform_vector(const VECTOR2& v) const
{
  VECTOR2 result = ROT2::invert(r) * v;

  #ifndef NEXCEPT
  boost::shared_ptr<const POSE2> pose;
  try
  {
    pose = shared_from_this();
  }
  catch (boost::bad_weak_ptr e)
  {
    std::cerr << "Pose2 can't check frames when allocated on the stack! Allocate on heap or disable debugging exceptions." << std::endl;
  }
  if (v.pose != rpose)
    throw FrameException();
  result.pose = boost::const_pointer_cast<const POSE2>(pose);
  #endif

  return result;
}

/// Transforms a point from one pose to another 
VECTOR2 POSE2::transform_point(const VECTOR2& p) const
{
  #ifndef NEXCEPT
  boost::shared_ptr<const POSE2> pose;
  try
  {
    pose = shared_from_this();
  }
  catch (boost::bad_weak_ptr e)
  {
    std::cerr << "Pose2 can't check frames when allocated on the stack! Allocate on heap or disable debugging exceptions." << std::endl;
  }
  if (p.pose != pose)
    throw FrameException();
  #endif

  VECTOR2 result = r * p + x;
  result.pose = rpose;
  return result;
}

/// Transforms a point from one pose to another 
VECTOR2 POSE2::inverse_transform_point(const VECTOR2& p) const
{
  VECTOR2 result = ROT2::invert(r) * (p - x);

  #ifndef NEXCEPT
  boost::shared_ptr<const POSE2> pose;
  try
  {
    pose = shared_from_this();
  }
  catch (boost::bad_weak_ptr e)
  {
    std::cerr << "Pose2 can't check frames when allocated on the stack! Allocate on heap or disable debugging exceptions." << std::endl;
  }
  if (p.pose != rpose)
    throw FrameException();
  result.pose = boost::const_pointer_cast<const POSE2>(pose);
  #endif

  return result;
}

/// Applies this pose to a vector 
VECTOR2 POSE2::transform_vector(boost::shared_ptr<const POSE2> target, const VECTOR2& v) 
{
  boost::shared_ptr<const POSE2> source = v.pose;

  // compute the relative transform
  TRANSFORM2 Tx = calc_transform(source, target);

  return Tx.transform_vector(v);
}

/// Transforms a point from one pose to another 
VECTOR2 POSE2::transform_point(boost::shared_ptr<const POSE2> target, const VECTOR2& point)
{
  boost::shared_ptr<const POSE2> source = point.pose;

  // compute the relative transform
  TRANSFORM2 Tx = calc_transform(source, target);

  // do the transform
  return Tx.transform_point(point);
}

/// Special method for inverting a 2D pose in place
POSE2& POSE2::invert()
{
  // get the new rotation 
  r.invert();

  // determine the new translation
  x = r * (-x);

  return *this;
}

/// Special method for inverting a 2D pose 
POSE2 POSE2::invert(const POSE2& p)
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
TRANSFORM2 POSE2::calc_transform(boost::shared_ptr<const POSE2> source, boost::shared_ptr<const POSE2> target)
{
  TRANSFORM2 result;
  boost::shared_ptr<const POSE2> r, s; 

  // setup the source and target
  result.source = source;
  result.target = target;

  // check for special case: no transformation 
  if (source == target)
  {
    result.r.set_identity();
    result.x.set_zero();
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
    result.r = left_r;      
    result.x = left_x;
    return result;
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

    // setup the pose 
    result.r = inv_right_r;      
    result.x = inv_right_x;
    return result;
  }

  // if both transforms are defined relative to the same frame, this is easy
  if (source->rpose == target->rpose)
  {
    // compute the inverse pose of p 
    ROT2 target_r = ROT2::invert(target->r);

    // multiply the inverse pose of p by this 
    result.r = target_r * source->r;      
    result.x = target_r * (source->x - target->x);
    return result;
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

    // multiply the inverse pose of p by this 
    result.r = inv_right_r * left_r;      
    result.x = inv_right_r * (left_x - right_x);
  }

  return result;
}

/// Outputs this matrix to the stream
std::ostream& Ravelin::operator<<(std::ostream& out, const POSE2& p)
{
  out << "r: " << p.r << " " << p.x << std::endl;
   
  return out;
}

