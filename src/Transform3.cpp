/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

/// Default constructor
/**
 * Sets matrix to the identity matrix
 */
TRANSFORM3::TRANSFORM3()
{
  set_identity();
}

/// Constructs a transformation from a unit quaternion and translation vector
TRANSFORM3::TRANSFORM3(const QUAT& q, const ORIGIN3& v)
{
  set(QUAT::normalize(q), v);
}

/// Constructs a transformation from a unit quaternion (for rotation) and zero translation
TRANSFORM3::TRANSFORM3(const QUAT& q)
{
  set(QUAT::normalize(q), ORIGIN3::zero());
}

/// Constructs a transformation from a rotation matrix and translation vector
TRANSFORM3::TRANSFORM3(const MATRIX3& r, const ORIGIN3& v) 
{
  set(r, v);
}

/// Constructs a transformation from a rotation matrix and zero translation
TRANSFORM3::TRANSFORM3(const MATRIX3& r)
{
  set(r, ORIGIN3::zero());
}

/// Constructs a transformation from a axis-angle representation and a translation vector
TRANSFORM3::TRANSFORM3(const AANGLE& a, const ORIGIN3& v)
{
  set(a, v);
}

/// Constructs a transformation from a axis-angle representation (for rotation) and zero translation
TRANSFORM3::TRANSFORM3(const AANGLE& a)
{
  set(a, ORIGIN3::zero());
}

/// Constructs a transformation using identity orientation and a translation vector
TRANSFORM3::TRANSFORM3(const ORIGIN3& v)
{
  set(QUAT::identity(), v);
}

TRANSFORM3& TRANSFORM3::operator=(const TRANSFORM3& p)
{
  q = p.q;
  x = p.x;
  source = p.source;
  target = p.target;

  return *this;
}

/// Determines whether two transformations in 3D are relatively equivalent
bool TRANSFORM3::rel_equal(const TRANSFORM3& p1, const TRANSFORM3& p2, REAL tol)
{
  // verify both sources and both targets are equal
  if (!(p1.source == p2.source && p1.target == p2.target))
    throw FrameException();

  // check x components first
  if (!OPS::rel_equal(p1.x[0], p2.x[0], tol) || 
      !OPS::rel_equal(p1.x[1], p2.x[1], tol) ||
      !OPS::rel_equal(p1.x[2], p2.x[2], tol))
    return false;

  // now check rotation (it's more expensive)
  REAL angle = QUAT::calc_angle(p1.q, p2.q);
  return OPS::rel_equal(angle, (REAL) 0.0, tol);
}

// Sets the transformation from an axis-angle representation and a translation vector
/**
 * \todo make this marginally faster by setting the matrix components directly (eliminates redundancy in setting non-rotation and non-translational components)
 */
TRANSFORM3& TRANSFORM3::set(const AANGLE& a, const ORIGIN3& v)
{
  this->q = a;
  this->x = v;
  return *this;
}

/// Sets the transformation from a rotation matrix and translation vector
/**
 * \todo make this marginally faster by setting the matrix components directly (eliminates redundancy in setting non-rotation and non-translational components)
 */
TRANSFORM3& TRANSFORM3::set(const MATRIX3& m, const ORIGIN3& v)
{
  q = m;
  x = v;
  return *this;
}

/// Sets the transformation from a unit quaternion and translation vector
TRANSFORM3& TRANSFORM3::set(const QUAT& q, const ORIGIN3& v)
{
  this->q = q;
  this->x = v;
  return *this;
}

/// Sets this transformation from a unit quaternion only (translation will be zero'd) 
TRANSFORM3& TRANSFORM3::set(const QUAT& q)
{
  this->q = q;
  this->x = ORIGIN3::zero();
  return *this;
}

/// Sets this transformation from a axis-angle object only (translation will be zero'd) 
TRANSFORM3& TRANSFORM3::set(const AANGLE& a)
{
  this->q = a;
  this->x = ORIGIN3::zero();
  return *this;
}

/// Sets this transformation from a 3x3 rotation matrix only (translation will be zero'd) 
TRANSFORM3& TRANSFORM3::set(const MATRIX3& m)
{
  this->q = m;
  this->x = ORIGIN3::zero();
  return *this;
}

/// Sets this matrix to identity
TRANSFORM3& TRANSFORM3::set_identity()
{
  this->q = QUAT::identity();
  this->x = ORIGIN3::zero();
  return *this;
}

/// Special method for inverting a transformation in place
TRANSFORM3& TRANSFORM3::invert()
{
  // get the new rotation 
  q.invert();

  // determine the new translation
  x = q * (-x);

  // swap the source and target poses
  std::swap(source, target);

  return *this;
}

/// Special method for inverting a transformation
TRANSFORM3 TRANSFORM3::invert(const TRANSFORM3& p)
{
  return TRANSFORM3(p).invert();
}

/// Multiplies this transformation by another
/**
 * If this transforms from frame b to frame c and T transforms from frame a 
 * to frame b, the multiplication transforms from from a to frame c.
 */
TRANSFORM3 TRANSFORM3::operator*(const TRANSFORM3& T) const
{
  #ifndef NEXCEPT
  if (source != T.target)
    throw FrameException();
  #endif
  TRANSFORM3 result(q * T.q, q*T.x + x);
  result.source = T.source;
  result.target = target;
  return result;
}

/// Transforms a pose
POSE3 TRANSFORM3::transform(const POSE3& p) const
{
  #ifndef NEXCEPT
  if (p.rpose != source)
    throw FrameException();
  #endif

  POSE3 result(target);
  result.q = q * p.q;
  result.x = ORIGIN3(q * p.x) + x; 
 
  return result;
}

/// Transforms a pose
POSE3 TRANSFORM3::inverse_transform(const POSE3& p) const
{
  #ifndef NEXCEPT
  if (p.rpose != target)
    throw FrameException();
  #endif

  POSE3 result(source);
  QUAT qi = QUAT::invert(q);
  result.q = qi * p.q;
  result.x = ORIGIN3(qi * p.x) - qi*x; 
 
  return result;
}

/// Transforms a vector from one pose to another 
VECTOR3 TRANSFORM3::transform(const VECTOR3& v) const
{
  #ifndef NEXCEPT
  if (v.pose != source)
    throw FrameException();
  #endif

  return VECTOR3(q * ORIGIN3(v), target);
}

/// Transforms a vector from one pose to another 
VECTOR3 TRANSFORM3::inverse_transform(const VECTOR3& v) const
{
  #ifndef NEXCEPT
  if (v.pose != target)
    throw FrameException();
  #endif

  return VECTOR3(QUAT::invert(q) * ORIGIN3(v), source);
}

/// Transforms a point from one pose to another 
POINT3 TRANSFORM3::transform(const POINT3& p) const
{
  #ifndef NEXCEPT
  if (p.pose != source)
    throw FrameException();
  #endif

  return POINT3(q * ORIGIN3(p) + x, target);
}

/// Transforms a point from one pose to another 
POINT3 TRANSFORM3::inverse_transform(const POINT3& p) const
{
  #ifndef NEXCEPT
  if (p.pose != target)
    throw FrameException();
  #endif

  return POINT3(QUAT::invert(q) * (ORIGIN3(p) - x), source);
}

/// Transforms a wrench from one pose to another 
WRENCH TRANSFORM3::transform(const WRENCH& w) const
{
  #ifndef NEXCEPT
  if (w.pose != source)
    throw FrameException();
  #endif

  // setup r and E
  ORIGIN3 r = -x;
  MATRIX3 E = q;
  VECTOR3 rv(r, target);
  const MATRIX3 ET = MATRIX3::transpose(E);

  // get the components of w
  ORIGIN3 top(w.get_force());
  ORIGIN3 bottom(w.get_torque());

  // do the calculations
  VECTOR3 Etop(E * top, target);
  VECTOR3 cross = VECTOR3::cross(rv, Etop);
  return WRENCH(Etop, (E * bottom) - cross, target);
}

/// Transforms a point from one pose to another 
WRENCH TRANSFORM3::inverse_transform(const WRENCH& w) const
{
  #ifndef NEXCEPT
  if (w.pose != target)
    throw FrameException();
  #endif

  // setup r and E
  MATRIX3 E = QUAT::invert(q);
  ORIGIN3 r = E * x;
  const MATRIX3 ET = MATRIX3::transpose(E);
  VECTOR3 rv(r, source);

  // get the components of w
  ORIGIN3 top(w.get_force());
  ORIGIN3 bottom(w.get_torque());

  // do the calculations
  VECTOR3 Etop(E * top, source);
  VECTOR3 cross = VECTOR3::cross(rv, Etop);
  return WRENCH(Etop, (E * bottom) - cross, source);
}

/// Transforms a twist from one pose to another 
TWIST TRANSFORM3::transform(const TWIST& t) const
{
  #ifndef NEXCEPT
  if (t.pose != source)
    throw FrameException();
  #endif

  // setup r and E
  ORIGIN3 r = -x;
  MATRIX3 E = q;
  VECTOR3 rv(r, target);
  const MATRIX3 ET = MATRIX3::transpose(E);

  // get the components of t 
  ORIGIN3 top(t.get_angular());
  ORIGIN3 bottom(t.get_linear());

  // do the calculations
  VECTOR3 Etop(E * top, target);
  VECTOR3 cross = VECTOR3::cross(rv, Etop);
  return TWIST(Etop, (E * bottom) - cross, target);
}

/// Transforms a twist from one pose to another 
TWIST TRANSFORM3::inverse_transform(const TWIST& t) const
{
  #ifndef NEXCEPT
  if (t.pose != target)
    throw FrameException();
  #endif

  // setup r and E
  MATRIX3 E = QUAT::invert(q);
  VECTOR3 r(E * x, source);
  const MATRIX3 ET = MATRIX3::transpose(E);

  // get the components of t 
  ORIGIN3 top(t.get_angular());
  ORIGIN3 bottom(t.get_linear());

  // do the calculations
  VECTOR3 Etop(E * top, source);
  VECTOR3 cross = VECTOR3::cross(r, Etop);
  return TWIST(Etop, (E * bottom) - cross, source);
}

/// Transforms a rigid body inertia from one pose to another 
SPATIAL_RB_INERTIA TRANSFORM3::transform(const SPATIAL_RB_INERTIA& J) const
{
  #ifndef NEXCEPT
  if (J.pose != source)
    throw FrameException();
  #endif

  // setup r and E
  ORIGIN3 r = -x;
  MATRIX3 E = q;
  VECTOR3 rv(r, target);
  const MATRIX3 ET = MATRIX3::transpose(E);

  // precompute some things
  VECTOR3 mr = rv * J.m;
  MATRIX3 rx = MATRIX3::skew_symmetric(rv);
  MATRIX3 hx = MATRIX3::skew_symmetric(J.h);
  MATRIX3 mrxrx = rx * MATRIX3::skew_symmetric(mr);  
  MATRIX3 EhxETrx = E * hx * ET * rx;

  // setup the new inertia
  SPATIAL_RB_INERTIA Jx(target);
  Jx.m = J.m;
  Jx.J = EhxETrx + MATRIX3::transpose(EhxETrx) + (E*J.J*ET) - mrxrx; 
  Jx.h = E * ORIGIN3(J.h) - mr;

  return Jx;
}

/// Transforms a rigid body inertia from one pose to another 
SPATIAL_RB_INERTIA TRANSFORM3::inverse_transform(const SPATIAL_RB_INERTIA& J) const
{
  #ifndef NEXCEPT
  if (J.pose != target)
    throw FrameException();
  #endif

  // setup r and E
  MATRIX3 E = QUAT::invert(q);
  ORIGIN3 r = E * x;
  VECTOR3 rv(r, target);
  const MATRIX3 ET = MATRIX3::transpose(E);

  // precompute some things
  VECTOR3 mr = rv * J.m;
  MATRIX3 rx = MATRIX3::skew_symmetric(rv);
  MATRIX3 hx = MATRIX3::skew_symmetric(J.h);
  MATRIX3 mrxrx = rx * MATRIX3::skew_symmetric(mr);  
  MATRIX3 EhxETrx = E * hx * ET * rx;

  // setup the new inertia
  SPATIAL_RB_INERTIA Jx(source);
  Jx.m = J.m;
  Jx.J = EhxETrx + MATRIX3::transpose(EhxETrx) + (E*J.J*ET) - mrxrx; 
  Jx.h = E * ORIGIN3(J.h) - mr;

  return Jx;
}

/// Transforms an articulated body inertia from one pose to another 
SPATIAL_AB_INERTIA TRANSFORM3::transform(const SPATIAL_AB_INERTIA& J) const
{
  #ifndef NEXCEPT
  if (J.pose != source)
    throw FrameException();
  #endif

  // setup r and E
  ORIGIN3 r = -x;
  MATRIX3 E = q;
  const MATRIX3 ET = MATRIX3::transpose(E);
  VECTOR3 rv(r, target);

  // precompute some things we'll need
  MATRIX3 rx = MATRIX3::skew_symmetric(rv);
  MATRIX3 HT = MATRIX3::transpose(J.H);
  MATRIX3 EJET = E * J.J * ET;
  MATRIX3 rx_E_HT_ET = rx*E*HT*ET;
  MATRIX3 EHET = E * J.H * ET;
  MATRIX3 EMET = E * J.M * ET;
  MATRIX3 rxEMET = rx * EMET;

  SPATIAL_AB_INERTIA result(target);
  result.M = EMET;
  result.H = EHET - rxEMET;
  result.J = EJET - rx_E_HT_ET + ((EHET - rxEMET) * rx); 

  return result;
}

/// Transforms an articulated body inertia from one pose to another 
SPATIAL_AB_INERTIA TRANSFORM3::inverse_transform(const SPATIAL_AB_INERTIA& J) const
{
  #ifndef NEXCEPT
  if (J.pose != target)
    throw FrameException();
  #endif

  // setup r and E
  MATRIX3 E = QUAT::invert(q);
  ORIGIN3 r = E * x;
  const MATRIX3 ET = MATRIX3::transpose(E);
  VECTOR3 rv(r, source);

  // precompute some things we'll need
  MATRIX3 rx = MATRIX3::skew_symmetric(rv);
  MATRIX3 HT = MATRIX3::transpose(J.H);
  MATRIX3 EJET = E * J.J * ET;
  MATRIX3 rx_E_HT_ET = rx*E*HT*ET;
  MATRIX3 EHET = E * J.H * ET;
  MATRIX3 EMET = E * J.M * ET;
  MATRIX3 rxEMET = rx * EMET;

  // compute the result
  SPATIAL_AB_INERTIA result(source);
  result.M = EMET;
  result.H = EHET - rxEMET;
  result.J = EJET - rx_E_HT_ET + ((EHET - rxEMET) * rx); 

  return result;
}

/// Outputs this matrix to the stream
std::ostream& Ravelin::operator<<(std::ostream& out, const TRANSFORM3& m)
{
  out << "q: " << m.q << " " << m.x << std::endl;
   
  return out;
}

