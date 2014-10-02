/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
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

/// Determines whether the transform is identity
bool TRANSFORM3::is_identity() const
{
  // first check the rotational component
  if (std::fabs(std::fabs(q.w) - (REAL) 1.0) > EPS)
    return false;

  // now check the translational component
  return x.norm() < EPS;
}

TRANSFORM3& TRANSFORM3::operator=(const TRANSFORM3& p)
{
  q = p.q;
  x = p.x;
  source = p.source;
  target = p.target;

  return *this;
}

/// Tranforms a vector with an interpolated transform (between transforms T1 and T2)
/**
 * \param T1 the pose to use when t=0
 * \param T2 the pose to use when t=1
 * \param t interpolation value
 * \param o the vector to transform 
 * \return the transformed vector 
 */
ORIGIN3 TRANSFORM3::interpolate_transform_vector(const TRANSFORM3& T1, const TRANSFORM3& T2, REAL t, const ORIGIN3& o)
{
  // interpolate the rotations
  QUAT q = QUAT::slerp(T1.q, T2.q, t);

  // transform the vector
  return q*o;
}

/// Tranforms a point with an interpolated pose (between poses T1 and T2)
/**
 * \param T1 the pose to use when t=0
 * \param T2 the pose to use when t=1
 * \param t interpolation value
 * \param o the point to transform 
 * \return the transformed point 
 */
ORIGIN3 TRANSFORM3::interpolate_transform_point(const TRANSFORM3& T1, const TRANSFORM3& T2, REAL t, const ORIGIN3& o)
{
  // interpolate the positions
  ORIGIN3 x = T1.x*(1-t) + T2.x*t;

  // interpolate the rotations
  QUAT q = QUAT::slerp(T1.q, T2.q, t);

  // trnsform the point 
  return x+q*o;
}

/// Tranforms a vector with the inverse of an interpolated transform (between transforms T1 and T2)
/**
 * \param T1 the pose to use when t=0
 * \param T2 the pose to use when t=1
 * \param t interpolation value
 * \param o the vector to transform 
 * \return the transformed vector 
 */
ORIGIN3 TRANSFORM3::interpolate_inverse_transform_vector(const TRANSFORM3& T1, const TRANSFORM3& T2, REAL t, const ORIGIN3& o)
{
  // interpolate the rotations
  QUAT q = QUAT::invert(QUAT::slerp(T1.q, T2.q, t));

  // transform the vector
  return q*o;
}

/// Tranforms a point with the inverse of an interpolated pose (between poses T1 and T2)
/**
 * \param T1 the pose to use when t=0
 * \param T2 the pose to use when t=1
 * \param t interpolation value
 * \param o the point to transform 
 * \return the transformed point 
 */
ORIGIN3 TRANSFORM3::interpolate_inverse_transform_point(const TRANSFORM3& T1, const TRANSFORM3& T2, REAL t, const ORIGIN3& o)
{
  // interpolate the rotations
  QUAT q = QUAT::invert(QUAT::slerp(T1.q, T2.q, t));

  // interpolate the positions
  ORIGIN3 x = q * (T1.x*(t-1) - T2.x*t);

  // trnsform the point 
  return x + q*o;
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
VECTOR3 TRANSFORM3::transform_vector(const VECTOR3& v) const
{
  #ifndef NEXCEPT
  if (v.pose != source)
    throw FrameException();
  #endif

  return VECTOR3(q * ORIGIN3(v), target);
}

/// Transforms a vector from one pose to another 
VECTOR3 TRANSFORM3::inverse_transform_vector(const VECTOR3& v) const
{
  #ifndef NEXCEPT
  if (v.pose != target)
    throw FrameException();
  #endif

  return VECTOR3(QUAT::invert(q) * ORIGIN3(v), source);
}

/// Transforms a point from one pose to another 
VECTOR3 TRANSFORM3::transform_point(const VECTOR3& p) const
{
  #ifndef NEXCEPT
  if (p.pose != source)
    throw FrameException();
  #endif

  return VECTOR3(q * ORIGIN3(p) + x, target);
}

/// Transforms a point from one pose to another 
VECTOR3 TRANSFORM3::inverse_transform_point(const VECTOR3& p) const
{
  #ifndef NEXCEPT
  if (p.pose != target)
    throw FrameException();
  #endif

  return VECTOR3(QUAT::invert(q) * (ORIGIN3(p) - x), source);
}

/// Transforms a force from one pose to another 
SFORCE TRANSFORM3::transform(const SFORCE& w) const
{
  SFORCE result;
  transform_spatial(w, result);
  return result;
}

/// Transforms a spatial vector
void TRANSFORM3::transform_spatial(const SVECTOR6& w, SVECTOR6& result) const
{
  #ifndef NEXCEPT
  if (w.pose != source)
    throw FrameException();
  #endif

  // setup r and E
  MATRIX3 E = q;
  ORIGIN3 r = E.transpose_mult(-x);
  VECTOR3 rv(r, source);

  // get the components of w
  VECTOR3 top = w.get_upper();
  VECTOR3 bottom = w.get_lower();

  // do the calculations
  VECTOR3 Etop(E * ORIGIN3(top), target);
  VECTOR3 cross = VECTOR3::cross(rv, top);
  result.set_upper(Etop);
  result.set_lower(VECTOR3(E * ORIGIN3(bottom - cross), target));
  result.pose = target;
}

/// Transforms a force from one pose to another 
SFORCE TRANSFORM3::inverse_transform(const SFORCE& w) const
{
  SFORCE result;
  inverse_transform_spatial(w, result);
  return result; 
}

/// Transforms a spatial vector from one pose to another 
void TRANSFORM3::inverse_transform_spatial(const SVECTOR6& w, SVECTOR6& result) const
{
  #ifndef NEXCEPT
  if (w.pose != target)
    throw FrameException();
  #endif

  // setup r and E
  MATRIX3 E = QUAT::invert(q);
  const ORIGIN3& r = x;
  VECTOR3 rv(r, source);

  // get the components of w
  VECTOR3 top = w.get_upper();
  VECTOR3 bottom = w.get_lower();

  // do the calculations
  VECTOR3 Etop(E * ORIGIN3(top), source);
  VECTOR3 cross = VECTOR3::cross(rv, top);
  result.set_upper(Etop);
  result.set_lower(VECTOR3(E * ORIGIN3(bottom - cross), source));
  result.pose = source;
}

/// Transforms a velocity from one pose to another 
SVELOCITY TRANSFORM3::transform(const SVELOCITY& t) const
{
  SVELOCITY result;
  transform_spatial(t, result);
  return result;
}

/// Transforms a velocity from one pose to another 
SVELOCITY TRANSFORM3::inverse_transform(const SVELOCITY& t) const
{
  SVELOCITY result;
  inverse_transform_spatial(t, result);
  return result;
}

/// Transforms a momentum from one pose to another 
SMOMENTUM TRANSFORM3::transform(const SMOMENTUM& t) const
{
  SMOMENTUM result;
  transform_spatial(t, result);
  return result;
}

/// Transforms a momentum from one pose to another 
SMOMENTUM TRANSFORM3::inverse_transform(const SMOMENTUM& t) const
{
  SMOMENTUM result;
  inverse_transform_spatial(t, result);
  return result;
}

/// Transforms an acceleration from one pose to another 
SACCEL TRANSFORM3::transform(const SACCEL& t) const
{
  SACCEL result;
  transform_spatial(t, result);
  return result;
}

/// Transforms an acceleration from one pose to another 
SACCEL TRANSFORM3::inverse_transform(const SACCEL& t) const
{
  SACCEL result;
  inverse_transform_spatial(t, result);
  return result;
}

/// Transforms a rigid body inertia from one pose to another 
/**
 * The operations for this come from:
 * | E     0 |  *  | -m*hx       eye(3)*m |  *  | E'    0  |        
 * | -E*rx E |     | J - m*hx*hx m*hx     |     | rx*E' E' |
 * which yields:
 * | A  m  |
 * | B  A' |
 * where:
 * A = -m*E*hx*E' + m*E*rx*E' 
 * B = m*E*rx*hx*E' + E*J*E' - m*E*hx*hx*E' - m*E*rx*rx*E' + m*E*hx*rx*E'
 */
SPATIAL_RB_INERTIA TRANSFORM3::transform(const SPATIAL_RB_INERTIA& J) const
{
  #ifndef NEXCEPT
  if (J.pose != source)
    throw FrameException();
  #endif

  // get r and E 
  MATRIX3 E = q;
  ORIGIN3 r = E.transpose_mult(-x);

  // precompute some things
  ORIGIN3 y = J.h - r;
//  MATRIX3 rx = MATRIX3::skew_symmetric(r);
//  MATRIX3 hx = MATRIX3::skew_symmetric(J.h);
//  MATRIX3 yx = MATRIX3::skew_symmetric(y);
//  MATRIX3 Z = J.J + (rx*hx) + (yx*rx);
  MATRIX3 Z = J.J;

  // transform the inertia 
  SPATIAL_RB_INERTIA Jx(target);
  Jx.m = J.m;
  Jx.h = E*y;
  Jx.J = E*Z*MATRIX3::transpose(E);

  return Jx;
}

/// Transforms a rigid body inertia from one pose to another 
SPATIAL_RB_INERTIA TRANSFORM3::inverse_transform(const SPATIAL_RB_INERTIA& J) const
{
  #ifndef NEXCEPT
  if (J.pose != target)
    throw FrameException();
  #endif

  // get r and E 
  MATRIX3 E = QUAT::invert(q);
  const ORIGIN3& r = x;

  // precompute some things
  ORIGIN3 y = J.h - r*J.m;
  MATRIX3 rx = MATRIX3::skew_symmetric(r);
  MATRIX3 hx = MATRIX3::skew_symmetric(J.h);
  MATRIX3 yx = MATRIX3::skew_symmetric(y);
  MATRIX3 Z = J.J + (rx*hx) + (yx*rx);

  // transform the inertia 
  SPATIAL_RB_INERTIA Jx(source);
  Jx.m = J.m;
  Jx.h = E*y;
  Jx.J = E*Z*MATRIX3::transpose(E);

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
  MATRIX3 E = q;
  const MATRIX3 ET = MATRIX3::transpose(E);
  ORIGIN3 r = ET.mult(-x);

  // precompute some things we'll need
  MATRIX3 rx = MATRIX3::skew_symmetric(r);
  MATRIX3 Y = J.H - rx*J.M;
  MATRIX3 EYET = E * Y * ET;
  MATRIX3 HT = MATRIX3::transpose(J.H);
  MATRIX3 Z = J.J - rx*HT + Y*rx;

  // setup the spatial inertia
  SPATIAL_AB_INERTIA result(target);
  result.M = E*J.M*ET;
  result.H = EYET;
  result.J = E*Z*ET; 

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
  const ORIGIN3& r = x;
  const MATRIX3 ET = MATRIX3::transpose(E);

  // precompute some things we'll need
  MATRIX3 rx = MATRIX3::skew_symmetric(r);
  MATRIX3 HT = MATRIX3::transpose(J.H);
  MATRIX3 EJET = E * J.J * ET;
  MATRIX3 EHET = E * J.H * ET;
  MATRIX3 EMET = E * J.M * ET;
  MATRIX3 rx_EMET = rx * EMET;
  MATRIX3 E_rx = E * rx;
  MATRIX3 E_rx_HT_ET = E_rx * HT * ET;
  MATRIX3 E_rx_M_rx_ET = E_rx * J.M * rx * ET;

  // setup the spatial inertia
  SPATIAL_AB_INERTIA result(source);
  result.M = EMET;
  result.H = EHET - rx_EMET;
  result.J = EJET - E_rx_HT_ET + EHET - E_rx_M_rx_ET; 

  return result;
}

/// Outputs this matrix to the stream
std::ostream& Ravelin::operator<<(std::ostream& out, const TRANSFORM3& m)
{
  out << "orientation: " << AANGLE(m.q) << " origin: " << m.x;
   
  return out;
}

