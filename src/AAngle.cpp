/// Default constructor
AANGLE::AANGLE()
{
  x = y = z = angle = 0;
}

/// Copy constructor
AANGLE::AANGLE(const QUAT& q)
{
  *this = q;
}

/// Constructs an axis-angle object from four REAL values
/*
 * \note automatically normalizes the axis
 */
AANGLE::AANGLE(REAL x, REAL y, REAL z, REAL angle)
{
  set(VECTOR3(x, y, z), angle);
}

/// Constructs an axis-angle object from a vector and a REAL value
/*
 * \note automatically normalizes the axis
 */
AANGLE::AANGLE(const VECTOR3& v, REAL angle)
{
  set(v, angle);
}

/// Constructs an axis-angle object from a 3x3 rotation matrix
AANGLE::AANGLE(const MATRIX3& m)
{
  *this = m;
}

/// Constructs an axis-angle object from a 3x3 rotation matrix and a desired axis
/*
 * \note automatically normalizes the axis
 */
AANGLE::AANGLE(const MATRIX3& m, const VECTOR3& v)
{
  set(m, v);
}

/// Sets this object from a three-dimensional vector and a REAL value
/**
 * \note automatically normalizes the axis
 */
AANGLE& AANGLE::set(const VECTOR3& v, REAL angle)
{
  VECTOR3 w = VECTOR3::normalize(v);
  x = w[0];
  y = w[1];
  z = w[2];
  this->angle = angle;
  return *this;
}

/// Sets the object from a rotation matrix and rotation axis
/**
 * \note automatically normalizes the axis
 */
void AANGLE::set(const MATRIX3& m, const VECTOR3& axis)
{
  // normalize the axis
  VECTOR3 axisn = VECTOR3::normalize(axis);

  // get the arc-cosine of the angle (clip it as necessary)
  REAL acosangle = (m.xx() + m.yy() + m.zz() - 1)/2; 
  if (acosangle > (REAL) 1.0)
    acosangle = (REAL) 1.0;
  else if (acosangle < -(REAL) 1.0)
    acosangle = -(REAL) 1.0;

  // compute angle of rotation
  angle = std::acos(acosangle);

  // set axis of rotation
  x = axisn.x();
  y = axisn.y();
  z = axisn.z();

  // if the angle of rotation is zero, can return now
  if (angle == (REAL) 0.0)
    return;

  // if the angle of rotation is PI, solve for axis
  REAL x,y,z;
  if (angle == M_PI)
  {
    x = safe_sqrt((m.xx()+1)/2);
    y = safe_sqrt((m.yy()+1)/2);
    z = safe_sqrt((m.zz()+1)/2);
  }
  else
  {
    REAL constant = (REAL) 1.0/((REAL) 2.0*std::sin(angle));
    x = constant * (m.zy() - m.yz());
    y = constant * (m.xz() - m.zx());
    z = constant * (m.yx() - m.xy());
  }

  // make sure that x, y, and z are not NaN
  assert(!std::isnan(x));
  assert(!std::isnan(y));
  assert(!std::isnan(z));

  // get the length of determined axis [x y z]; if length is zero, angle is
  // zero...
  REAL len = safe_sqrt(x*x + y*y + z*z);
  if (len == (REAL) 0.0)
  {
    angle = (REAL) 0.0;
    return;
  }

  // normalize vector [x y z] (generally not necessary, but safe...)
  if (std::fabs(len - (REAL) 1.0) > EPS)
  {
    REAL ilen = (REAL) 1.0/len;
    x *= ilen;
    y *= ilen;
    z *= ilen;
  }

  // form the determined axis and verify it points in same or exact opposite
  // direction as axis
  VECTOR3 determine_axis(x,y,z, axis.pose);
  
  // reverse the angle, if necessary
  if (VECTOR3::dot(determine_axis, axisn) < 0)
    angle = -angle;
}

/// Copies a axis-angle object
AANGLE& AANGLE::operator=(const AANGLE& a)
{
  x = a.x;
  y = a.y;
  z = a.z;
  angle = a.angle;
  return *this;
}

/// Multiplies two axis-angle representations
/**
 * This has the effect of performing a rotation by a, then performing a 
 * rotation by <b>this</b>.
 */
AANGLE AANGLE::operator*(const AANGLE& a) const
{
  QUAT q1 = *this;
  QUAT q2 = a;
  q1 *= q2;
  return AANGLE(q1);
}

/// Multiplies this axis-angle representation by another and stores the result in <b>this</b>
/**
 * This has the effect of performing a rotation by a, then performing a
 * rotation by <b>this</b>.
 */
AANGLE& AANGLE::operator*=(const AANGLE& a)
{
  QUAT q1 = *this;
  QUAT q2 = a;
  q1 *= q2;
  *this = q1;
  return *this;
}

/// Sends the representation to the specified stream
std::ostream& Ravelin::operator<<(std::ostream& out, const AANGLE& a)
{
  out << "[ " << a.x << ' ' << a.y << ' ' << a.z << " ] "  << a.angle << "  ";
  return out;
}

/// Sets axis angle from unit quaternion 
/**
 * \todo test this method
 */
AANGLE& AANGLE::operator=(const QUAT& q)
{
  QUAT qn = QUAT::normalize(q);
  AANGLE& a = *this;
  a.x = qn.x;
  a.y = qn.y;
  a.z = qn.z;
  a.angle = std::acos(qn.w) * 2.0f;

  // verify that [x y z] normalized
  REAL nrm = AANGLE::safe_sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
  if (std::fabs(nrm) < EPS)
  {
    // axis is zero; set it to [1 0 0] arbitrarily
    a.x = 1.0f;
    a.y = a.z = 0.0f;
  }
  else if (std::fabs(nrm - 1.0f) > EPS)
  {
    a.x /= nrm;
    a.y /= nrm;
    a.z /= nrm;
  }

  return a;
}

/// Sets the object from a rotation matrix
AANGLE& AANGLE::operator=(const MATRIX3& m)
{
  // get the arc-cosine of the angle (clip it as necessary)
  REAL acosangle = (m.xx() + m.yy() + m.zz() - 1)* (REAL) 0.5; 
  if (acosangle > (REAL) 1.0)
    acosangle = (REAL) 1.0;
  else if (acosangle < (REAL) -1.0)
    acosangle = (REAL) -1.0;

  // compute angle of rotation
  AANGLE& a = *this;
  a.angle = std::acos(acosangle);
  
  // if angle is 0, then axis is arbitrary
  if (a.angle < std::numeric_limits<REAL>::epsilon())
  {
    a.x = (REAL) 1.0;
    a.y = a.z = (REAL) 0.0;
    return a;
  }
  
  // if angle is pi then must solve for rx, ry, rz
  if (std::fabs(a.angle-M_PI) < std::numeric_limits<REAL>::epsilon())
  {
    a.x = AANGLE::safe_sqrt((m.xx()+1)*(REAL) 0.5);
    a.y = AANGLE::safe_sqrt((m.yy()+1)*(REAL) 0.5);
    a.z = AANGLE::safe_sqrt((m.zz()+1)*(REAL) 0.5);
    assert(!std::isnan(a.x));
    assert(!std::isnan(a.y));
    assert(!std::isnan(a.z));
    return a;
  }
  
  // standard case
  REAL constant = (REAL) 1.0/((REAL) 2.0*std::sin(a.angle));
  a.x = constant * (m.zy() - m.yx());
  a.y = constant * (m.xz() - m.zx());
  a.z = constant * (m.yx() - m.xy());
  
  // normalize the axis (generally not necessary, but safe...)
  REAL len = AANGLE::safe_sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
  assert(len != (REAL) 0.0);
  if (std::fabs(len - (REAL) 1.0) > EPS)
  {
    REAL ilen = (REAL) 1.0/len;
    a.x *= ilen;
    a.y *= ilen;
    a.z *= ilen;
  }

  assert(!std::isnan(a.angle));
  assert(!std::isnan(a.x));
  assert(!std::isnan(a.y));
  assert(!std::isnan(a.z));

  return a;
}

