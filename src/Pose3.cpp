/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

using boost::shared_ptr;

/// Default constructor
/**
 * Sets matrix to the identity matrix
 */
POSE3::POSE3(boost::shared_ptr<const POSE3> relative_pose)
{
  set_identity();
  rpose = relative_pose;
}

/// Constructs a pose from a unit quaternion and translation vector
POSE3::POSE3(const QUAT& q, const ORIGIN3& v, boost::shared_ptr<const POSE3> relative_pose)
{
  set(QUAT::normalize(q), v);
  rpose = relative_pose;
}

/// Constructs a pose from a unit quaternion (for rotation) and zero translation
POSE3::POSE3(const QUAT& q, boost::shared_ptr<const POSE3> relative_pose)
{
  set(QUAT::normalize(q), ORIGIN3::zero());
  rpose = relative_pose;
}

/// Constructs a pose from a rotation matrix and translation vector
POSE3::POSE3(const MATRIX3& r, const ORIGIN3& v, boost::shared_ptr<const POSE3> relative_pose) 
{
  set(r, v);
  rpose = relative_pose;
}

/// Constructs a pose from a rotation matrix and zero translation
POSE3::POSE3(const MATRIX3& r, boost::shared_ptr<const POSE3> relative_pose)
{
  set(r, ORIGIN3::zero());
  rpose = relative_pose;
}

/// Constructs a pose from a axis-angle representation and a translation vector
POSE3::POSE3(const AANGLE& a, const ORIGIN3& v, boost::shared_ptr<const POSE3> relative_pose)
{
  set(a, v);
  rpose = relative_pose;
}

/// Constructs a pose from a axis-angle representation (for rotation) and zero translation
POSE3::POSE3(const AANGLE& a, boost::shared_ptr<const POSE3> relative_pose)
{
  set(a, ORIGIN3::zero());
  rpose = relative_pose;
}

/// Constructs a pose using identity orientation and a translation vector
POSE3::POSE3(const ORIGIN3& v, boost::shared_ptr<const POSE3> relative_pose)
{
  set(QUAT::identity(), v);
  rpose = relative_pose;
}

POSE3& POSE3::operator=(const POSE3& p)
{
  q = p.q;
  x = p.x;
  rpose = p.rpose;
  return *this;
}

/// Determines whether two poses in 3D are relatively equivalent
bool POSE3::rel_equal(const POSE3& p1, const POSE3& p2, REAL tol)
{
  #ifndef NEXCEPT
  // verify relative frames are identical 
  if (p1.rpose != p2.rpose)
    throw FrameException();
  #endif

  // check x components first
  if (!OPS::rel_equal(p1.x[0], p2.x[0], tol) || 
      !OPS::rel_equal(p1.x[1], p2.x[1], tol) ||
      !OPS::rel_equal(p1.x[2], p2.x[2], tol))
    return false;

  // now check rotation (it's more expensive)
  REAL angle = QUAT::calc_angle(p1.q, p2.q);
  return OPS::rel_equal(angle, (REAL) 0.0, tol);
}

/// Tranforms a vector with an interpolated pose (between poses P1 and P2)
/**
 * \param P1 the pose to use when t=0
 * \param P2 the pose to use when t=1
 * \param t interpolation value
 * \param o the vector to transform 
 * \return the transformed vector 
 */
VECTOR3 POSE3::interpolate_transform_vector(const POSE3& P1, const POSE3& P2, REAL t, const ORIGIN3& o)
{
  #ifndef NEXCEPT
  if (P1.rpose != P2.rpose)
    throw FrameException();
  #endif

  // interpolate the rotations
  QUAT q = QUAT::slerp(P1.q, P2.q, t);

  // transform the vector
  return VECTOR3(q*o, P1.rpose);
}

/// Tranforms a point with an interpolated pose (between poses P1 and P2)
/**
 * \param P1 the pose to use when t=0
 * \param P2 the pose to use when t=1
 * \param t interpolation value
 * \param o the point to transform 
 * \return the transformed point 
 */
VECTOR3 POSE3::interpolate_transform_point(const POSE3& P1, const POSE3& P2, REAL t, const ORIGIN3& o)
{
  #ifndef NEXCEPT
  if (P1.rpose != P2.rpose)
    throw FrameException();
  #endif

  // interpolate the positions
  ORIGIN3 x = P1.x*(1-t) + P2.x*t;

  // interpolate the rotations
  QUAT q = QUAT::slerp(P1.q, P2.q, t);

  // trnsform the point 
  return VECTOR3(x+q*o, P1.rpose);
}

/// Multiplies the quaternion derivative
VECTOR3 POSE3::qG_mult(REAL qdw, REAL qdx, REAL qdy, REAL qdz) const
{
  const shared_ptr<const POSE3> GLOBAL;

  // this will put everything into the global frame
  TRANSFORM3 T = calc_relative_pose(shared_from_this(), GLOBAL);

  // use the resulting q
  return T.q.G_mult(qdx, qdy, qdz, qdw);
}

/// Multiplies the quaternion derivative
QUAT POSE3::qG_transpose_mult(const VECTOR3& omega) const
{
  const shared_ptr<const POSE3> GLOBAL;

  // this will put everything into the global frame
  TRANSFORM3 T = calc_relative_pose(shared_from_this(), GLOBAL);

  // use the resulting q
  return T.q.G_transpose_mult(omega);
}

// Sets the pose from an axis-angle representation and a translation vector
/**
 * \note relative pose will be unaltered
 */
POSE3& POSE3::set(const AANGLE& a, const ORIGIN3& v)
{
  this->q = a;
  this->x = v;
  return *this;
}

/// Sets the pose from a rotation matrix and translation vector
/**
 * \note relative pose will be unaltered
 */
POSE3& POSE3::set(const MATRIX3& m, const ORIGIN3& v)
{
  q = m;
  x = v;
  return *this;
}

/// Sets the pose from a unit quaternion and translation vector
/**
 * \note relative pose will be unaltered
 */
POSE3& POSE3::set(const QUAT& q, const ORIGIN3& v)
{
  this->q = q;
  this->x = v;
  return *this;
}

/// Sets this pose from a unit quaternion only (translation will be zero'd) 
/**
 * \note relative pose will be unaltered
 */
POSE3& POSE3::set(const QUAT& q)
{
  this->q = q;
  this->x = ORIGIN3::zero();
  return *this;
}

/// Sets this pose from a axis-angle object only (translation will be zero'd) 
/**
 * \note relative pose will be unaltered
 */
POSE3& POSE3::set(const AANGLE& a)
{
  this->q = a;
  this->x = ORIGIN3::zero();
  return *this;
}

/// Sets this pose from a 3x3 rotation matrix only (translation will be zero'd) 
/**
 * \note relative pose will be unaltered
 */
POSE3& POSE3::set(const MATRIX3& m)
{
  this->q = m;
  this->x = ORIGIN3::zero();
  return *this;
}

/// Sets this matrix to identity
/**
 * \note relative pose will be unaltered
 */
POSE3& POSE3::set_identity()
{
  this->q = QUAT::identity();
  this->x = ORIGIN3::zero();
  return *this;
}

/// method for inverting a pose in place
POSE3& POSE3::invert()
{
  // get the new rotation 
  q.invert();

  // determine the new translation
  x = q * (-x);

  return *this;
}

/// Special method for inverting a pose
POSE3 POSE3::invert(const POSE3& p)
{
  return POSE3(p).invert();
}

/// Multiplies this pose by another
POSE3 POSE3::operator*(const POSE3& p) const
{
  #ifndef NEXCEPT
  boost::shared_ptr<const POSE3> pose;
  try
  {
    pose = shared_from_this();
  }
  catch (boost::bad_weak_ptr e)
  {
    std::cerr << "Pose3::operator*() - pose allocated on stack!" << std::endl;
  }
  if (p.rpose != pose)
    throw FrameException();
  #endif

  return POSE3(q * p.q, q*p.x + x, rpose);
}

/// Transforms a vector from one pose to another 
VECTOR3 POSE3::inverse_transform_vector(const VECTOR3& v) const
{
  boost::shared_ptr<const POSE3> pose;
  try
  {
    pose = shared_from_this();
  }
  catch (boost::bad_weak_ptr e)
  {
    std::cerr << "Pose3::inverse_transform() - pose allocated on stack!" << std::endl;
  }
  #ifndef NEXCEPT
  if (v.pose != rpose)
    throw FrameException();
  #endif

  return transform_vector(pose, v);
}

/// Transforms a point from one pose to another 
VECTOR3 POSE3::inverse_transform_point(const VECTOR3& p) const
{
  boost::shared_ptr<const POSE3> pose;
  try
  {
    pose = shared_from_this();
  }
  catch (boost::bad_weak_ptr e)
  {
    std::cerr << "Pose3::inverse_transform() - pose allocated on stack!" << std::endl;
  }
  #ifndef NEXCEPT
  if (p.pose != rpose)
    throw FrameException();
  #endif

  return transform_point(pose, p); 
}

/// Transforms a point from one pose to another 
SFORCE POSE3::inverse_transform(const SFORCE& w) const
{
  // get the target pose
  boost::shared_ptr<const POSE3> pose;
  try
  {
    pose = shared_from_this();
  }
  catch (boost::bad_weak_ptr e)
  {
    std::cerr << "Pose3::inverse_transform() - pose allocated on stack!" << std::endl;
  }
  #ifndef NEXCEPT
  if (w.pose != rpose)
    throw FrameException();
  #endif

  return transform(pose, w);
}

/// Transforms a vector of forcees 
std::vector<SFORCE>& POSE3::transform(boost::shared_ptr<const POSE3> target, const std::vector<SFORCE>& w, std::vector<SFORCE>& result)
{
  // look for empty vector (easy case)
  if (w.empty())
  {
    result.clear();
    return result;
  }

  // setup the source pose
  boost::shared_ptr<const POSE3> source = w[0].pose; 

  #ifndef NEXCEPT
  for (unsigned i=1; i< w.size(); i++)
    if (source != w[i].pose)
      throw FrameException();
  #endif

  // quick check
  if (source == target)
    return (result = w);

  // compute the relative transform
  TRANSFORM3 Tx = calc_transform(source, target);

  // setup r and E
  VECTOR3 r;
  MATRIX3 E;
  get_r_E(Tx, r, E);

  // resize the result vector
  result.resize(w.size());

  // look over all forcees
  for (unsigned i=0; i< w.size(); i++)
    transform_spatial(target, w[i], r, E, result[i]);

  return result;
}

/// transforms a spatial vector using precomputation
void POSE3::transform_spatial(boost::shared_ptr<const POSE3> target, const SVECTOR6& w, const VECTOR3& r, const MATRIX3& E, SVECTOR6& result)
{
  // get the components of w[i] 
  VECTOR3 top = w.get_upper();
  VECTOR3 bottom = w.get_lower();

  // do the calculations
  VECTOR3 Etop(E * ORIGIN3(top), target);
  VECTOR3 cross = VECTOR3::cross(r, top);
  result.set_upper(Etop);
  result.set_lower(VECTOR3(E * ORIGIN3(bottom - cross), target));
  result.pose = target;
} 

/// Transforms the force 
SFORCE POSE3::transform(boost::shared_ptr<const POSE3> target, const SFORCE& v)
{
  // setup the force
  SFORCE f;

  // do the transform
  transform_spatial(target, v, f);

  return f;
}

/// Transforms a spatial vector
void POSE3::transform_spatial(boost::shared_ptr<const POSE3> target, const SVECTOR6& v, SVECTOR6& s)
{
  // setup the source pose 
  boost::shared_ptr<const POSE3> source = v.pose;

  // quick check
  if (source == target)
  {
    s=v;
    return;
  }

  // compute the relative transform
  TRANSFORM3 Tx = calc_transform(source, target);

  // setup r and E
  VECTOR3 r;
  MATRIX3 E;
  get_r_E(Tx, r, E);

  // get the components of v
  VECTOR3 top = v.get_upper();
  VECTOR3 bottom = v.get_lower();

  // do the calculations
  VECTOR3 Etop(E * ORIGIN3(top), target);
  VECTOR3 cross = VECTOR3::cross(r, top);
  s.set_upper(Etop);
  s.set_lower(VECTOR3(E * ORIGIN3(bottom - cross), target));
  s.pose = target;
}

/// Transforms an acceleration from one pose to another 
SACCEL POSE3::inverse_transform(const SACCEL& t) const
{
  // determine the resulting frame
  boost::shared_ptr<const POSE3> pose;
  try
  {
    pose = shared_from_this();
  }
  catch (boost::bad_weak_ptr e)
  {
    std::cerr << "Pose3::inverse_transform() - pose allocated on stack!" << std::endl;
  }
  #ifndef NEXCEPT
  if (t.pose != rpose)
    throw FrameException();
  #endif

  return transform(pose, t); 
}

/// Transforms the acceleration 
SACCEL POSE3::transform(boost::shared_ptr<const POSE3> target, const SACCEL& t)
{
  // setup the acceleration 
  SACCEL a;

  // do the transform
  transform_spatial(target, t, a);

  return a;
}

/// Transforms a vector of accelerations 
std::vector<SACCEL>& POSE3::transform(boost::shared_ptr<const POSE3> target, const std::vector<SACCEL>& t, std::vector<SACCEL>& result)
{
  // look for empty vector (easy case)
  if (t.empty())
  {
    result.clear();
    return result;
  }

  // setup the source pose
  boost::shared_ptr<const POSE3> source = t[0].pose; 

  #ifndef NEXCEPT
  for (unsigned i=1; i< t.size(); i++)
    if (source != t[i].pose)
      throw FrameException();
  #endif

  // quick check
  if (source == target)
    return (result = t);

  // compute the relative transform
  TRANSFORM3 Tx = calc_transform(source, target);

  // setup r and E
  VECTOR3 r;
  MATRIX3 E;
  get_r_E(Tx, r, E);

  // the spatial transformation is:
  // | E    0 |
  // | Erx' E |

  // resize the result vector
  result.resize(t.size());

  // transform 
  for (unsigned i=0; i< t.size(); i++)
    transform_spatial(target, t[i], r, E, result[i]);

  return result;
}

/// Transforms a velocity from one pose to another 
SVELOCITY POSE3::inverse_transform(const SVELOCITY& t) const
{
  // determine the resulting frame
  boost::shared_ptr<const POSE3> pose;
  try
  {
    pose = shared_from_this();
  }
  catch (boost::bad_weak_ptr e)
  {
    std::cerr << "Pose3::inverse_transform() - pose allocated on stack!" << std::endl;
  }
  #ifndef NEXCEPT
  if (t.pose != rpose)
    throw FrameException();
  #endif

  return transform(pose, t);
}

/// Transforms the velocity  
SVELOCITY POSE3::transform(boost::shared_ptr<const POSE3> target, const SVELOCITY& t)
{
  SVELOCITY v;
  transform_spatial(target, t, v);
  return v;
}

/// Transforms a vector of velocities 
std::vector<SVELOCITY>& POSE3::transform(boost::shared_ptr<const POSE3> target, const std::vector<SVELOCITY>& t, std::vector<SVELOCITY>& result)
{
  // look for empty vector (easy case)
  if (t.empty())
  {
    result.clear();
    return result;
  }

  // setup the source pose
  boost::shared_ptr<const POSE3> source = t[0].pose; 

  #ifndef NEXCEPT
  for (unsigned i=1; i< t.size(); i++)
    if (source != t[i].pose)
      throw FrameException();
  #endif

  // quick check
  if (source == target)
    return (result = t);

  // compute the relative transform
  TRANSFORM3 Tx = calc_transform(source, target);

  // setup r and E
  VECTOR3 r;
  MATRIX3 E;
  get_r_E(Tx, r, E);

  // resize the result vector
  result.resize(t.size());

  // transform the individual vectors 
  for (unsigned i=0; i< t.size(); i++)
    transform_spatial(target, t[i], r, E, result[i]);

  return result;
}

/// Transforms a momentum from one pose to another 
SMOMENTUM POSE3::inverse_transform(const SMOMENTUM& t) const
{
  // determine the resulting frame
  boost::shared_ptr<const POSE3> pose;
  try
  {
    pose = shared_from_this();
  }
  catch (boost::bad_weak_ptr e)
  {
    std::cerr << "Pose3::inverse_transform() - pose allocated on stack!" << std::endl;
  }
  #ifndef NEXCEPT
  if (t.pose != rpose)
    throw FrameException();
  #endif

  return transform(pose, t); 
}

/// Transforms the momentum 
SMOMENTUM POSE3::transform(boost::shared_ptr<const POSE3> target, const SMOMENTUM& t)
{
  SMOMENTUM result;
  transform_spatial(target, t, result);
  return result;
}

/// Transforms a vector of momenta 
std::vector<SMOMENTUM>& POSE3::transform(boost::shared_ptr<const POSE3> target, const std::vector<SMOMENTUM>& t, std::vector<SMOMENTUM>& result)
{
  // look for empty vector (easy case)
  if (t.empty())
  {
    result.clear();
    return result;
  }

  // setup the source pose
  boost::shared_ptr<const POSE3> source = t[0].pose; 

  // quick check
  if (source == target)
    return (result = t);

  #ifndef NEXCEPT
  for (unsigned i=1; i< t.size(); i++)
    if (source != t[i].pose)
      throw FrameException();
  #endif

  // compute the relative transform
  TRANSFORM3 Tx = calc_transform(source, target);

  // setup r and E
  VECTOR3 r;
  MATRIX3 E;
  get_r_E(Tx, r, E);

  // resize the result vector
  result.resize(t.size());

  // transform all momenta 
  for (unsigned i=0; i< t.size(); i++)
    transform_spatial(target, t[i], r, E, result[i]);

  return result;
}

/// Transforms a rigid body inertia from one pose to another 
SPATIAL_RB_INERTIA POSE3::inverse_transform(const SPATIAL_RB_INERTIA& J) const
{
  boost::shared_ptr<const POSE3> pose;
  try
  {
    pose = shared_from_this();
  }
  catch (boost::bad_weak_ptr e)
  {
    std::cerr << "Pose3::inverse_transform() - pose allocated on stack!" << std::endl;
  }
  #ifndef NEXCEPT
  if (J.pose != rpose)
    throw FrameException();
  #endif

  return transform(pose, J);
}

/// Transforms an articulated body inertia from one pose to another 
SPATIAL_AB_INERTIA POSE3::inverse_transform(const SPATIAL_AB_INERTIA& J) const
{
  boost::shared_ptr<const POSE3> pose;
  try
  {
    pose = shared_from_this();
  }
  catch (boost::bad_weak_ptr e)
  {
    std::cerr << "Pose3::inverse_transform() - pose allocated on stack!" << std::endl;
  }
  #ifndef NEXCEPT
  if (J.pose != rpose)
    throw FrameException();
  #endif

  return transform(pose, J);
}

/// Gets r and E from a transform 
/**
 * \param T a transformation from frame A to frame B
 * \param r on return, the vector from A's origin to B's origin in A frame 
 * \param E rotates vectors in A's orientation to vectors in B's orientation 
 */ 
void POSE3::get_r_E(const TRANSFORM3& T, VECTOR3& r, MATRIX3& E)
{
  // x is the translation from frame A to frame B

  // p_a
  // bTa * p_a = bQa * p_a + bXa

  // note that x is translation from relative pose to this pose
  // q is rotation from vectors in this pose to relative pose 
  E = T.q;
  r = VECTOR3(E.transpose_mult(-T.x), T.source);
}

/// Gets r and E from the current pose only
/**
 * Assume the current pose is frame A and the relative pose is frame B.
 * \param r on return, the vector from A's origin to B's origin in A frame 
 * \param E rotates vectors in A's orientation to vectors in B's orientation 
 */ 
void POSE3::get_r_E(VECTOR3& r, MATRIX3& E, bool inverse) const
{
  // note that x is translation from relative pose to this pose
  // q is rotation from vectors in this pose to relative pose 
  if (!inverse)
  {
    E = q;
    r = VECTOR3(E.transpose_mult(-x), shared_from_this());
  }
  else
  {
    E = QUAT::invert(q);
    r = VECTOR3(x, rpose);
  } 
}

/// Applies this pose to a vector 
VECTOR3 POSE3::transform_vector(boost::shared_ptr<const POSE3> target, const VECTOR3& v) 
{
  // setup source pose 
  boost::shared_ptr<const POSE3> source = v.pose;

  // quick check
  if (source == target)
    return v;

  // do the transformation
  return calc_transform(source, target).transform_vector(v);
}

/// Transforms a point from one pose to another 
VECTOR3 POSE3::transform_point(boost::shared_ptr<const POSE3> target, const VECTOR3& point)
{
  // setup source pose 
  boost::shared_ptr<const POSE3> source = point.pose;

  #ifndef NEXCEPT
  if (source != point.pose)
    throw FrameException();
  #endif

  // quick check
  if (source == target)
    return point;

  // do the transformation
  return calc_transform(source, target).transform_point(point);
}

SMOMENTUM POSE3::transform(const SMOMENTUM& t) const { return transform(rpose, t); }
SFORCE POSE3::transform(const SFORCE& w) const { return transform(rpose, w); }
SVELOCITY POSE3::transform(const SVELOCITY& t) const { return transform(rpose, t); }
SACCEL POSE3::transform(const SACCEL& t) const { return transform(rpose, t); }

/// Transforms a spatial articulated body inertia 
SPATIAL_AB_INERTIA POSE3::transform(boost::shared_ptr<const POSE3> target, const SPATIAL_AB_INERTIA& m)
{
  // setup source pose 
  boost::shared_ptr<const POSE3> source = m.pose;

  // quick check
  if (source == target)
    return m;

  // compute the transformation
  return calc_transform(source, target).transform(m);
}

/// Transforms a spatial RB inertia to the given pose
SPATIAL_RB_INERTIA POSE3::transform(boost::shared_ptr<const POSE3> target, const SPATIAL_RB_INERTIA& J)
{
  // setup source pose 
  boost::shared_ptr<const POSE3> source = J.pose;

  // quick check
  if (source == target)
    return J;

  // do the transformation
  return calc_transform(source, target).transform(J);
}

/// Determines whether pose p exists in the chain of relative poses (and, if so, how far down the chain)
bool POSE3::is_common(boost::shared_ptr<const POSE3> x, boost::shared_ptr<const POSE3> p, unsigned& i)
{
  // reset i
  i = 0;

  while (true)
  {
    if (x == p)
      return true;
    if (!x)
      return false;
    x = x->rpose;
    i++;
  }
}

/// Computes the relative transformation from this pose to another
TRANSFORM3 POSE3::calc_transform(boost::shared_ptr<const POSE3> source, boost::shared_ptr<const POSE3> target)
{
  boost::shared_ptr<const POSE3> r, s; 
  TRANSFORM3 result;

  // setup the source and targets
  result.source = source;
  result.target = target;

  // check for special case: no transformation 
  if (source == target)
  {
    result.q.set_identity();
    result.x.set_zero();
    return result;
  }

  // check for special case: transformation to global frame
  if (!target)
  {
    // combine transforms from this to i: this will give aTl
    result.q = source->q;
    result.x = source->x;
    s = source;
    while (s)
    {
      s = s->rpose;
      if (!s)
        break;
      result.x = s->x + s->q * result.x;
      result.q = s->q * result.q;
    }

    // setup the transform 
    return result; 
  }

  // check for special case: transformation from global frame
  if (!source)
  {
    // combine transforms from target to q
    result.q = target->q;
    result.x = target->x;
    r = target;
    while (r)
    {
      r = r->rpose;
      if (!r)
        break;
      result.x = r->x + r->q * result.x; 
      result.q = r->q * result.q;
    }

    // compute the inverse pose of the right 
    result.q = QUAT::invert(result.q);
    result.x = result.q * (-result.x);
    return result;
  }

  // if both transforms are defined relative to the same frame, this is easy
  if (source->rpose == target->rpose)
  {
    // compute the inverse pose of p 
    QUAT inv_target_q = QUAT::invert(target->q);
    result.q = inv_target_q * source->q;
    result.x = inv_target_q * (source->x - target->x);
    return result;
  }
  else
  {
    // search for the common link - we arbitrary move up the target while
    // one step at a time while searching through all levels of the source 
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
    
    // combine transforms from this to i: this will give rTs, where r is
    // the common frame and s is the source
    QUAT left_q = source->q;
    ORIGIN3 left_x = source->x;
    s = source;
    for (unsigned j=0; j < i; j++)
    {
      s = s->rpose;
      if (!s)
        break;
      left_x = s->x + s->q * left_x;
      left_q = s->q * left_q;
    }

    // combine transforms from target to r
    QUAT right_q = target->q;
    ORIGIN3 right_x = target->x;
    while (target != r)
    {
      target = target->rpose;
      if (!target)
        break;
      right_x = target->x + target->q * right_x; 
      right_q = target->q * right_q;
    }

    // compute the inverse pose of the right 
    QUAT inv_right_q = QUAT::invert(right_q);

    // multiply the inverse pose of p by this 
    result.q = inv_right_q * left_q;      
    result.x = inv_right_q * (left_x - right_x);
    return result;
  }
}

/// Adds a velocity to a pose to yield a new pose
POSE3 POSE3::operator+(const SVELOCITY& v) const
{
  #ifndef NEXCEPT
  if (rpose != v.pose)
    throw FrameException();
  #endif

  // setup the new pose
  POSE3 P;
  P.rpose = rpose;
  P.x = x + ORIGIN3(v.get_linear());
  P.q = q + QUAT::deriv(q, v.get_angular());
  P.q.normalize();
  return P;
}

/// Computes the spatial velocity that gets from one pose (P1) to another (P2)
SVELOCITY POSE3::diff(const POSE3& P1, const POSE3& P2)
{
  #ifndef NEXCEPT
  if (P1.rpose != P2.rpose)
    throw FrameException();
  #endif

  // compute the differential in positions
  VECTOR3 xd(P2.x - P1.x, P1.rpose);

  // compute the differential in orientations
  VECTOR3 omega = QUAT::to_omega(P1.q, P2.q - P1.q);
  omega.pose = P1.rpose;

  return SVELOCITY(xd, omega, P1.rpose);
}

/// Updates the relative pose for this pose while retaining the same absolute pose
POSE3& POSE3::update_relative_pose(boost::shared_ptr<const POSE3> pose)
{
  // compute the transform from the current relative pose to the new relative
  // pose
  TRANSFORM3 T = calc_transform(rpose, pose);

  // transform the quaternion
  q = T.q * q;

  // transform the origin 
  x = T.q * x + T.x;

  // update the relative pose
  rpose = pose;

  return *this;
}

/// Outputs this matrix to the stream
std::ostream& Ravelin::operator<<(std::ostream& out, const POSE3& m)
{
  out << "orientation: " << AANGLE(m.q) << " origin: " << m.x << " relative pose: " << m.rpose;
   
  return out;
}

