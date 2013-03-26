/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

/// Default constructor
/**
 * Sets matrix to the identity matrix
 */
POSE::POSE()
{
  set_identity();
}

/// Constructs a 4x4 homogeneous transformation matrix from a unit quaternion and translation vector
POSE::POSE(const QUAT& q, const VECTOR3& v)
{
  set(QUAT::normalize(q), v);
}

/// Constructs a 4x4 homogeneous transformation matrix from a unit quaternion (for rotation) and zero translation
POSE::POSE(const QUAT& q)
{
  set(QUAT::normalize(q), VECTOR3::zero());
}

/// Constructs a 4x4 homogeneous transformation matrix from a rotation matrix and translation vector
POSE::POSE(const MATRIX3& r, const VECTOR3& v) 
{
  set(r, v);
}

/// Constructs a 4x4 homogeneous transformation matrix from a rotation matrix and zero translation
POSE::POSE(const MATRIX3& r)
{
  set(r, VECTOR3::zero());
}

/// Constructs a 4x4 homogeneous transformation matrix from a axis-angle representation and a translation vector
POSE::POSE(const AANGLE& a, const VECTOR3& v)
{
  set(a, v);
}

/// Constructs a 4x4 homogeneous transformation matrix from a axis-angle representation (for rotation) and zero translation
POSE::POSE(const AANGLE& a)
{
  set(a, VECTOR3::zero());
}

/// Constructs a 4x4 homogeneous transformation matrix using identity orientation and a translation vector
POSE::POSE(const VECTOR3& v)
{
  set(QUAT::identity(), v);
}

POSE& POSE::operator=(const POSE& p)
{
  q = p.q;
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
POSE POSE::interpolate(const POSE& T1, const POSE& T2, REAL t)
{
  // interpolate the positions
  VECTOR3 x = T1.x*(1-t) + T2.x*t;

  // interpolate the rotations
  QUAT q = QUAT::slerp(T1.q, T2.q, t);

  // return the new matrix
  return POSE(q, x);
}

// Sets the matrix to be a 4x4 homogeneous transform from an axis-angle representation and a translation vector
/**
 * \todo make this marginally faster by setting the matrix components directly (eliminates redundancy in setting non-rotation and non-translational components)
 */
POSE& POSE::set(const AANGLE& a, const VECTOR3& v)
{
  this->q = a;
  this->x = v;
  return *this;
}

/// Sets the matrix to be a 4x4 homogeneous transform from a rotation matrix and translation vector
/**
 * \todo make this marginally faster by setting the matrix components directly (eliminates redundancy in setting non-rotation and non-translational components)
 */
POSE& POSE::set(const MATRIX3& m, const VECTOR3& v)
{
  q = m;
  x = v;
  return *this;
}

/// Sets the pose from a unit quaternion and translation vector
POSE& POSE::set(const QUAT& q, const VECTOR3& v)
{
  this->q = q;
  this->x = v;
  return *this;
}

/// Sets this pose from a unit quaternion only (translation will be zero'd) 
POSE& POSE::set(const QUAT& q)
{
  this->q = q;
  this->x = VECTOR3::zero();
  return *this;
}

/// Sets this pose from a axis-angle object only (translation will be zero'd) 
POSE& POSE::set(const AANGLE& a)
{
  this->q = a;
  this->x = VECTOR3::zero();
  return *this;
}

/// Sets this pose from a 3x3 rotation matrix only (translation will be zero'd) 
POSE& POSE::set(const MATRIX3& m)
{
  this->q = m;
  this->x = VECTOR3::zero();
  return *this;
}

/// Sets this matrix to identity
POSE& POSE::set_identity()
{
  this->q = QUAT::identity();
  this->x = VECTOR3::zero();
  return *this;
}

/// Applies this pose to a vector 
VECTOR3 POSE::mult_vector(const VECTOR3& v) const
{
  return q*v;
}

/// Applies the inverse of this pose to a vector 
VECTOR3 POSE::inverse_mult_vector(const VECTOR3& v) const
{
  return QUAT::invert(q) * v;
}

/// Transforms a point 
VECTOR3 POSE::mult_point(const VECTOR3& v) const
{
  return q*v + x;
}

/// Applies the inverse of this pose to a point 
VECTOR3 POSE::inverse_mult_point(const VECTOR3& v) const
{
  QUAT iq = QUAT::invert(q);
  return iq*v - iq*x;
}

/// Special method for inverting a 4x4 transformation matrix in place
POSE& POSE::invert()
{
  // get the new rotation 
  q.inverse();

  // determine the new translation
  x = q * (-x);

  return *this;
}

/// Special method for inverseing a 4x4 transformation matrix
POSE POSE::inverse(const POSE& p)
{
  return POSE(p).invert();
}

/// Multiplies this pose by another
POSE POSE::operator*(const POSE& p) const
{
  #ifndef NEXCEPT
  if (rpose != p.rpose)
    throw FrameException();
  #endif
  return POSE(q * p.q, q*p.x + x);
}

/// Transforms a spatial articulated body inertia 
SPATIAL_AB_INERTIA POSE::transform(boost::shared_ptr<const POSE> p, const SPATIAL_AB_INERTIA& m) const
{
  return transform(shared_from_this(), p, m);
}

/// Transforms a spatial articulated body inertia 
SPATIAL_AB_INERTIA POSE::transform(boost::shared_ptr<const POSE> source, boost::shared_ptr<const POSE> target, const SPATIAL_AB_INERTIA& m)
{
  #ifndef NEXCEPT
  if (source != m.pose)
    throw FrameException();
  #endif

  // compute the relative transform
  std::pair<QUAT, VECTOR3> Tx = calc_transform(source, target);

  // setup r and E
  const VECTOR3& r = Tx.second;
  const MATRIX3 E = Tx.first;
  const MATRIX3 ET = MATRIX3::transpose(E);

  // precompute some things we'll need
  MATRIX3 rx = MATRIX3::skew_symmetric(r);
  MATRIX3 HT = MATRIX3::transpose(m.H);
  MATRIX3 EJET = E * m.J * ET;
  MATRIX3 rx_E_HT_ET = rx*E*HT*ET;
  MATRIX3 EHET = E * m.H * ET;
  MATRIX3 EMET = E * m.M * ET;
  MATRIX3 rxEMET = rx * EMET;

  SPATIAL_AB_INERTIA result;
  result.pose = boost::const_pointer_cast<POSE>(target);
  result.M = EMET;
  result.H = EHET - rxEMET;
  result.J = EJET - rx_E_HT_ET + ((EHET - rxEMET) * rx); 
  return result;
}

/// Transforms the wrench 
WRENCH POSE::transform(boost::shared_ptr<const POSE> p, const WRENCH& v) const
{
  return transform(shared_from_this(), p, v);
}

/// Transforms the wrench 
WRENCH POSE::transform(boost::shared_ptr<const POSE> source, boost::shared_ptr<const POSE> target, const WRENCH& v)
{
  #ifndef NEXCEPT
  if (source != v.pose)
    throw FrameException();
  #endif

  // compute the relative transform
  std::pair<QUAT, VECTOR3> Tx = calc_transform(source, target);

  // setup r and E
  const VECTOR3& r = Tx.second;
  const MATRIX3 E = Tx.first;

  // get the components of v
  VECTOR3 top = v.get_force();
  VECTOR3 bottom = v.get_torque();

  // do the calculations
  VECTOR3 Etop = E * top;
  VECTOR3 cross = VECTOR3::cross(r, Etop);
  WRENCH result(Etop, (E * bottom) - cross);
  result.pose = boost::const_pointer_cast<POSE>(target);
  return result;
}

/// Transforms the twist 
TWIST POSE::transform(boost::shared_ptr<const POSE> p, const TWIST& t) const
{
  return transform(shared_from_this(), p, t);
}

/// Transforms the twist 
TWIST POSE::transform(boost::shared_ptr<const POSE> source, boost::shared_ptr<const POSE> target, const TWIST& t)
{
  #ifndef NEXCEPT
  if (source != t.pose)
    throw FrameException();
  #endif

  // compute the relative transform
  std::pair<QUAT, VECTOR3> Tx = calc_transform(source, target);

  // setup r and E
  const VECTOR3& r = Tx.second;
  const MATRIX3 E = Tx.first;

  // get the components of t 
  VECTOR3 top = t.get_angular();
  VECTOR3 bottom = t.get_linear();

  // do the calculations
  VECTOR3 Etop = E * top;
  VECTOR3 cross = VECTOR3::cross(r, Etop);
  TWIST result(Etop, (E * bottom) - cross);
  result.pose = boost::const_pointer_cast<POSE>(target);
  return result;
}

/// Transforms a spatial RB inertia to the given pose
SPATIAL_RB_INERTIA POSE::transform(boost::shared_ptr<const POSE> p, const SPATIAL_RB_INERTIA& J) const
{
  return transform(shared_from_this(), p, J);
}

/// Transforms a spatial RB inertia to the given pose
SPATIAL_RB_INERTIA POSE::transform(boost::shared_ptr<const POSE> source, boost::shared_ptr<const POSE> target, const SPATIAL_RB_INERTIA& J)
{
  #ifndef NEXCEPT
  if (source != J.pose)
    throw FrameException();
  #endif

  // compute the relative transform
  std::pair<QUAT, VECTOR3> Tx = calc_transform(source, target);

  // setup r and E
  const VECTOR3& r = Tx.second;
  const MATRIX3 E = Tx.first;
  const MATRIX3 ET = MATRIX3::transpose(E);

  // precompute some things
  VECTOR3 mr = r * J.m;
  MATRIX3 rx = MATRIX3::skew_symmetric(r);
  MATRIX3 hx = MATRIX3::skew_symmetric(J.h);
  MATRIX3 mrxrx = rx * MATRIX3::skew_symmetric(mr);  
  MATRIX3 EhxETrx = E * hx * ET * rx;

  // setup the new inertia
  SPATIAL_RB_INERTIA Jx;
  Jx.pose = boost::const_pointer_cast<POSE>(target);
  Jx.m = J.m;
  Jx.J = EhxETrx + MATRIX3::transpose(EhxETrx) + (E*J.J*ET) - mrxrx; 
  Jx.h = E * J.h - mr;
  return Jx;
}

/// Determines whether pose p exists in the chain of relative poses (and, if so, how far down the chain)
bool POSE::is_common(boost::shared_ptr<const POSE> x, boost::shared_ptr<const POSE> p, unsigned& i)
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
std::pair<QUAT, VECTOR3> POSE::calc_transform(boost::shared_ptr<const POSE> p) const
{
  return calc_transform(shared_from_this(), p);
}

/// Computes the relative transformation from this pose to another
std::pair<QUAT, VECTOR3> POSE::calc_transform(boost::shared_ptr<const POSE> source, boost::shared_ptr<const POSE> target)
{
  std::pair<QUAT, VECTOR3> result;
  boost::shared_ptr<const POSE> r, s; 

  // check for special case: transformation to and from global frame
  if (source == target && !source)
  {
    result.first = QUAT::identity();
    result.second.set_zero();
    return result;
  }

  // check for special case: transformation to global frame
  if (!target)
  {
    // combine transforms from this to i: this will give aTl
    QUAT left_q = source->q;
    VECTOR3 left_x = source->x;
    s = source;
    while (s)
    {
      s = s->rpose;
      if (!s)
        break;
      left_x = s->x + s->q * left_x;
      left_q = s->q * left_q;
    }

    // setup the transform 
    result.first = left_q;      
    result.second = -left_x;
  }

  // check for special case: transformation from global frame
  if (!source)
  {
    // combine transforms from target to q
    QUAT right_q = target->q;
    VECTOR3 right_x = target->x;
    r = target;
    while (r)
    {
      r = r->rpose;
      if (!r)
        break;
      right_x = r->x + r->q * right_x; 
      right_q = r->q * right_q;
    }

    // compute the inverse pose of the right 
    QUAT inv_right_q = QUAT::invert(right_q);
    VECTOR3 inv_right_x = inv_right_q * (-right_x);

    // multiply the inverse pose of p by this 
    result.first = inv_right_q;      
    result.second = inv_right_q * inv_right_x;
    return result;
  }

  // if both transforms are defined relative to the same frame, this is easy
  if (source->rpose == target->rpose)
  {
    // compute the inverse pose of p 
    QUAT target_q = QUAT::invert(target->q);
    VECTOR3 target_x = target_q * (-target->x);

    // multiply the inverse pose of p by this 
    result.first = target_q * source->q;      
    result.second = target_q * (target_x - source->x);
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
    QUAT left_q = source->q;
    VECTOR3 left_x = source->x;
    s = source;
    for (unsigned j=0; j < i; j++)
    {
      s = s->rpose;
      left_x = s->x + s->q * left_x;
      left_q = s->q * left_q;
    }

    // combine transforms from target to q
    QUAT right_q = target->q;
    VECTOR3 right_x = target->x;
    while (target != r)
    {
      target = target->rpose;
      right_x = target->x + target->q * right_x; 
      right_q = target->q * right_q;
    }

    // compute the inverse pose of the right 
    QUAT inv_right_q = QUAT::invert(right_q);
    VECTOR3 inv_right_x = inv_right_q * (-right_x);

    // multiply the inverse pose of p by this 
    result.first = inv_right_q * left_q;      
    result.second = inv_right_q * (inv_right_x - left_x);
  }

  return result;
}

/// Outputs this matrix to the stream
std::ostream& Ravelin::operator<<(std::ostream& out, const POSE& m)
{
  out << "q: " << m.q << " " << m.x << std::endl;
   
  return out;
}

