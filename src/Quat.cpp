/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

/// Default constructor
/**
 * Quaternion constructed as [0,0,0] 1
 */
QUAT::QUAT()
{
  x = y = z = (REAL) 0.0;
  w = (REAL) 1.0;
}

/// Constructs a quaternion from four REAL values
/**
 * \param x the first component of the vector
 * \param y the second component of the vector
 * \param z the third component of the vector
 * \param w the scalar component of the quaternion
 */
QUAT::QUAT(REAL x, REAL y, REAL z, REAL w)
{
  this->x = x;
  this->y = y;
  this->z = z;
  this->w = w;
}

/// Constructs a zero quaternion
QUAT QUAT::zero()
{
  QUAT q;
  q.x = q.y = q.z = q.w = (REAL) 0.0;
  return q;
}

/// Gets the i'th component of the quaternion
/**
 * \param i for i=0, returns the first component of the vector (x), for
 *        i=1, returns the second component of the vector (y), for i=2,
 *        returns the third component of the vector (z), for i=3, returns
 *        the scalar (w).
 */
REAL& QUAT::operator[](unsigned i)
{
  switch (i)
  {
    case 0: return x;
    case 1: return y;
    case 2: return z;
    case 3: return w;
    default:
      throw InvalidIndexException();
  }

  // make compiler happy
  return x;
}

/// Gets the i'th component of the quaternion
/**
 * \param i for i=0, returns the first component of the vector (x), for
 *        i=1, returns the second component of the vector (y), for i=2,
 *        returns the third component of the vector (z), for i=3, returns
 *        the scalar (w).
 */
const REAL& QUAT::operator[](unsigned i) const
{
  switch (i)
  {
    case 0: return x;
    case 1: return y;
    case 2: return z;
    case 3: return w;
    default:
      throw InvalidIndexException();
  }

  // make compiler happy
  return x;
}

/// Copy operator
QUAT& QUAT::operator=(const QUAT& q)
{
  x = q.x;
  y = q.y;
  z = q.z;
  w = q.w;
  return *this;
}

/// Negates the value of each of the x, y, and z coordinates in place
void QUAT::conjugate()
{
  *this = conjugate(*this);
}

/// Negates the value of each of the x, y, and z coordinates of the given quaternion
/**
 * \param q the quaternion from which the conjugate is to be computed
 */
QUAT QUAT::conjugate(const QUAT& q)
{
  return QUAT(-q.x, -q.y, -q.z, q.w);
}

/// Multiplies the 3x4 matrix 'L' times a four dimensional array
/**
 * This matrix is used in the relationships omega' = 2*L*qd and 
 * alpha' = 2*L*qdd, where omega'/alpha' are the angular velocity/acceleration 
 * of a rigid body in the body's frame and qd/qdd are the first/second time 
 * derivatives of the Euler (unit quaternion) parameters.
 */
VECTOR3 QUAT::L_mult(REAL qx, REAL qy, REAL qz, REAL qw) const
{
  const double e0 = qw;
  const double e1 = qx;
  const double e2 = qy;
  const double e3 = qz;

  VECTOR3 v;
  v.x() = -x*e0 + w*e1 + z*e2 - y*e3;
  v.y() = -y*e0 - z*e1 + w*e2 + x*e3;
  v.z() = -z*e0 + y*e1 - x*e2 + w*e3;
  return v;
}

/*
/// Computes the matrix 'L' used for generalized coordinate calculations
MatrixNf& QUAT::determine_L(MatrixNf& L) const
{
  L.resize(3,4);
  L(0,0) = -x;  L(0,1) = +w;  L(0,2) = +z;  L(0,3) = -y;
  L(1,0) = -y;  L(1,1) = -z;  L(1,2) = +w;  L(1,3) = +x;
  L(2,0) = -z;  L(2,1) = +y;  L(2,2) = -x;  L(2,3) = +w;
  return L;
}

/// Computes the matrix 'G' used for generalized coordinate calculations
MatrixNf& QUAT::determine_G(MatrixNf& G) const
{
  const unsigned X = 0, Y = 1, Z = 2, W = 3;
  G.resize(3,4);
  G(X,X) = -x;  G(X,Y) = +w;  G(X,Z) = -z;  G(X,W) = +y;
  G(Y,X) = -y;  G(Y,Y) = +z;  G(Y,Z) = +w;  G(Y,W) = -x;
  G(Z,X) = -z;  G(Z,Y) = -y;  G(Z,Z) = +x;  G(Z,W) = +w;
  return G;
}
*/

/// Multiplies the matrix 'G' by a quaternion vector
/**
 * This matrix is used in the relationships omega = 2*G*qd and 
 * alpha = 2*G*qdd, where omega/alpha are the angular velocity/acceleration 
 * of a rigid body in the game frame and qd/qdd are the first/second time 
 * derivatives of the Euler (unit quaternion) parameters.
 */
VECTOR3 QUAT::G_mult(REAL qx, REAL qy, REAL qz, REAL qw) const
{
  const double e0 = qw;
  const double e1 = qx;
  const double e2 = qy;
  const double e3 = qz;

  VECTOR3 r;
  r.x() = -x*e0 + w*e1 - z*e2 + y*e3;
  r.y() = -y*e0 + z*e1 + w*e2 - x*e3;
  r.z() = -z*e0 - y*e1 + x*e2 + w*e3;
  return r;
}

/// Multiplies the transpose of the matrix 'G' by a vector
/**
 * This matrix is used in the relationships qd = 1/2*G^T*omega and 
 * qdd = 1/2*G^T*alpha - 1/4*omega^2*q, where omega/alpha are the angular 
 * velocity/acceleration of a rigid body in the global frame and qd/qdd are the 
 * first/second time derivatives of the Euler (unit quaternion) parameters.
 */
QUAT QUAT::G_transpose_mult(const VECTOR3& v) const
{
  QUAT q;
  q.w = -x*v.x() - y*v.y() - z*v.z();
  q.x = +w*v.x() + z*v.y() - y*v.z();
  q.y = -z*v.x() + w*v.y() + x*v.z();
  q.z = +y*v.x() - x*v.y() + w*v.z();
  return q;
}

/// Multiplies the transpose of the matrix 'L' by a vector
/**
 * This matrix is used in the relationships qd = 1/2*L^T*omega and 
 * qdd = 1/2*L^T*alpha' - 1/4*omega'^2*q, where omega'/alpha' are the angular 
 * velocity/acceleration of a rigid body in the body frame and qd/qdd are the 
 * first/second time derivatives of the Euler (unit quaternion) parameters.
 */
QUAT QUAT::L_transpose_mult(const VECTOR3& v) const
{
  QUAT q;
  q.w = -x*v.x() - y*v.y() - z*v.z();
  q.x = +w*v.x() - z*v.y() + y*v.z();
  q.y = +z*v.x() + w*v.y() - x*v.z();
  q.z = -y*v.x() + x*v.y() + w*v.z();
  return q;
}

/// Computes the derivative of a quaternion given current orientation and angular velocity
/**
 * Note that the derivative is not generally a unit quaternion.
 * \param q the current orientation
 * \param w the angular velocity (in the global frame)
 * Uses the matrix: 
 *      |  -q.x  +q.w  -q.z  +q.y  |
 * G =  |  -q.y  +q.z  +q.w  -q.x  |
 *      |  -q.z  -q.y  +q.x  +q.w  |
 */
QUAT QUAT::deriv(const QUAT& q, const VECTOR3& w)
{
  QUAT qd;

  qd.w = .5 * (-q.x * w.x() - q.y * w.y() - q.z * w.z()); 
  qd.x = .5 * (+q.w * w.x() + q.z * w.y() - q.y * w.z());
  qd.y = .5 * (-q.z * w.x() + q.w * w.y() + q.x * w.z());
  qd.z = .5 * (+q.y * w.x() - q.x * w.y() + q.w * w.z());

  return qd;
}

/// Computes the angular velocity of a body given the current quaternion orientation and the quaternion velocity
VECTOR3 QUAT::to_omega(const QUAT& q, const QUAT& qd)
{
  VECTOR3 omega;
  omega.x() = 2 * (-q.x * qd.w + q.w * qd.x - q.z * qd.y + q.y * qd.z);
  omega.y() = 2 * (-q.y * qd.w + q.z * qd.x + q.w * qd.y - q.x * qd.z);
  omega.z() = 2 * (-q.z * qd.w - q.y * qd.x + q.x * qd.y + q.w * qd.z);
  return omega;
}

/// Calculates the second derivative of a quaternion
/**
 * \note alpha and omega are acceleration and velocity vectors in the global 
 *       frame
 */
QUAT QUAT::dderiv(const QUAT& q, const VECTOR3& omega, const VECTOR3& alpha)
{
  QUAT qdd = QUAT::deriv(q, alpha) - (REAL) 0.25 * omega.norm_sq() * q;
 
  return qdd;
}

/// Performs spherical linear interpolation between this and q
/**
 * Sets the orientation represented by this quaternion to a linear interpolated
 * orientation between this and q.  If alpha is 0, then *this is unaltered.
 * If alpha is 1, then *this is set to q.
 * \param q the "final" orientation
 * \param alpha interpolation value (0 <= alpha <= 1)
 */
void QUAT::slerp(const QUAT& q, REAL alpha)
{
  *this = slerp(*this, q, alpha);
}

/// Performs linear interpolation between this and q
/**
 * Sets the orientation represented by this quaternion to a linear interpolated
 * orientation between this and q.  If alpha is 0, then *this is unaltered.
 * If alpha is 1, then *this is set to q.
 * \param q the "final" orientation
 * \param alpha interpolation value (0 <= alpha <= 1)
 */
void QUAT::lerp(const QUAT& q, REAL alpha)
{
  *this = lerp(*this, q, alpha);
}

/// Performs non-spherical linear interpolation between this and q
/**
 * Calculates the orientation of determined by linear interpolation between
 * q1 and q2.  If alpha is 0, then q1 is returned.  If alpha is 1, then 
 * q2 is returned. 
 * \param q1 the "initial" orientation
 * \param q2 the "final" orientation
 * \param alpha interpolation value (0 <= alpha <= 1)
 * \return the orientation linearly interpolated between q1 and q2 
 */
QUAT QUAT::lerp(const QUAT& q1, const QUAT& q2, REAL alpha)
{
  if (alpha < (REAL) 0.0 || alpha > (REAL) 1.0)
    throw std::runtime_error("Attempting to interpolate using QUAT::lerp() with t not in interval [0,1]");

  // compute q1'q2
  REAL dot = q1.x*q2.x + q1.y*q2.y + q1.z*q2.z + q1.w*q2.w;

  // see whether we need to use the conjugate of q2
  if (dot < (REAL) 0.0)  
  {
    QUAT q = q1*((REAL) 1.0 - alpha) - q2*alpha;
    q.normalize();
    return q;
  }
  else
  {
    QUAT q = q1*((REAL) 1.0 - alpha) + q2*alpha;
    q.normalize();
    return q;
  }
}

/// Determines whether two poses in 3D are relatively equivalent
bool QUAT::rel_equal(const QUAT& q1, const QUAT& q2, REAL tol)
{
  // update the tolerance if necessary
  if (tol < (REAL) 0.0)
    tol = EPS;

  // now check rotation (it's more expensive)
  REAL angle = QUAT::calc_angle(q1, q2);
  return OPS::rel_equal(angle, (REAL) 0.0, tol);
}

/// Calculates the angle between two orientations
REAL QUAT::calc_angle(const QUAT& q1, const QUAT& q2)
{
  // compute q1'q2
  REAL dot = q1.x*q2.x + q1.y*q2.y + q1.z*q2.z + q1.w*q2.w;

  // clip dot
  if (dot < (REAL) -1.0)
    dot = (REAL) -1.0;
  else if (dot > (REAL) 1.0)
    dot = (REAL) 1.0;

  // find the angle between the two
  return std::acos(std::fabs(dot));
}

/// Sets up the quaternion from roll-pitch-yaw 
/**
 * \param roll the rotation around the local x axis
 * \param pitch the rotation around the local y axis
 * \param yaw the rotation around the local z axis
 */
QUAT QUAT::rpy(REAL roll, REAL pitch, REAL yaw)
{
  const REAL PHI = roll * (REAL) 0.5;
  const REAL THE = pitch * (REAL) 0.5;
  const REAL PSI = yaw * (REAL) 0.5;

  // precompute trig fns
  const REAL CPHI = std::cos(PHI);
  const REAL SPHI = std::sin(PHI);
  const REAL CPSI = std::cos(PSI);
  const REAL SPSI = std::sin(PSI);
  const REAL CTHE = std::cos(THE);
  const REAL STHE = std::sin(THE);

  // construct Quaternion
  QUAT q;
  q.w = CPHI * CTHE * CPSI + SPHI * STHE * SPSI;
  q.x = SPHI * CTHE * CPSI - CPHI * STHE * SPSI;
  q.y = CPHI * STHE * CPSI + SPHI * CTHE * SPSI;
  q.z = CPHI * CTHE * SPSI - SPHI * STHE * CPSI;
  return q;
}

/// Converts quaternion to roll/pitch/yaw 
void QUAT::to_rpy(REAL& roll, REAL& pitch, REAL& yaw) const
{
  pitch = std::atan2((REAL) 2.0*(y*z + w*x), w*w - x*x - y*y + z*z);
  yaw = std::asin((REAL) -2.0*(x*z - w*y));
  roll = std::atan2((REAL) 2.0*(x*y + w*z), w*w + x*x - y*y - z*z);
}

/// Performs spherical linear interpolation between this and q
/**
 * Calculates the orientation of determined by linear interpolation between
 * q1 and q2.  If alpha is 0, then q1 is returned.  If alpha is 1, then 
 * q2 is returned. 
 * \param q1 the "initial" orientation
 * \param q2 the "final" orientation
 * \param alpha interpolation value (0 <= alpha <= 1)
 * \return the orientation linearly interpolated between q1 and q2
 * \todo rewrite this function to avoid cancellation errors 
 */
QUAT QUAT::slerp(const QUAT& q1, const QUAT& q2, REAL alpha)
{
  if (alpha < (REAL) 0.0 || alpha > (REAL) 1.0)
    throw std::runtime_error("Attempting to interpolate using QUAT::slerp() with t not in interval [0,1]");

  // compute q1'q2
  REAL dot = q1.x*q2.x + q1.y*q2.y + q1.z*q2.z + q1.w*q2.w;

  // see whether we need to use the conjugate of q2
  bool use_conj = (dot < (REAL) 0.0);  

  // clip dot
  if (dot < (REAL) -1.0)
    dot = (REAL) -1.0;
  else if (dot > (REAL) 1.0)
    dot = (REAL) 1.0;

  // find the angle between the two
  REAL theta = std::acos(std::fabs(dot));
  if (theta == 0.0)
    return q1;

  // do slerp
  REAL sin1at = std::sin((1-alpha)*theta);
  REAL sinat = std::sin(alpha*theta);
  REAL sint_i = (REAL) 1.0/std::sin(theta);
  QUAT qa(q1.x*sin1at, q1.y*sin1at, q1.z*sin1at, q1.w*sin1at);
  QUAT qb(q2.x*sinat, q2.y*sinat, q2.z*sinat, q2.w*sinat);
  if (use_conj)
    qb.conjugate();
  QUAT qc = qa + qb;
  qc.x *= sint_i;
  qc.y *= sint_i;
  qc.z *= sint_i;
  qc.w *= sint_i;

  return qc; 
}

/// Computes the inverse orientation represented by this quaternion in place
QUAT& QUAT::invert()
{
  conjugate();
  return *this;
}

/// Computes the inverse orientation of a quaternion
QUAT QUAT::invert(const QUAT& q)
{
  return conjugate(q);
}

/// Normalizes this quaternion in place
void QUAT::normalize()
{
  REAL mag = magnitude();
  if (mag != (REAL) 0.0)
  {
    REAL magi = (REAL) 1.0/mag;
    x *= magi;
    y *= magi;
    z *= magi;
    w *= magi;
  }
  else
  {
    x = y = z = 0;
    w = 1;
  }
}

/// Computes the normalized quaternion of q
QUAT QUAT::normalize(const QUAT& q)
{
  QUAT qn(q);
  qn.normalize();
  return qn;
}

/// Determines whether this is a unit quaternion
bool QUAT::unit() const
{
  REAL mag = magnitude();
  return (std::fabs(mag - (REAL) 1.0) < std::sqrt(EPS));
}

/// Subtracts a quaternion from <b>this</b>
QUAT QUAT::operator-(const QUAT& q) const
{
  return QUAT(x+q.x, y+q.y, z+q.z, w+q.w);
}

/// Subtracts one quaternion from another 
QUAT& QUAT::operator-=(const QUAT& q)
{
  *this = *this - q;
  return *this;
}

/// Adds this quaternion to another and returns the result in a new quaternion
QUAT QUAT::operator+(const QUAT& q) const
{
  return QUAT(x+q.x, y+q.y, z+q.z, w+q.w);
}

/// Adds this quaon to another and stores the result in <b>this</b>
QUAT& QUAT::operator+=(const QUAT& q)
{
  *this = *this + q;
  return *this;
}

/// Multiplies inv(q) by <b>this</b> and returns the result
QUAT QUAT::operator/(const QUAT& q) const
{
  return conjugate(q) * (*this);
}

/// Multiplies <b>this</b> by q and returns the result 
QUAT QUAT::operator*(const QUAT& q) const
{
  QUAT qm;

  qm.w = w * q.w - x * q.x - y * q.y - z * q.z; 
  qm.x = w * q.x + x * q.w + y * q.z - z * q.y;
  qm.y = w * q.y + y * q.w + z * q.x - x * q.z;
  qm.z = w * q.z + z * q.w + x * q.y - y * q.x;

  return qm;
}

/// Multiplies a quaternion by a scalar
QUAT QUAT::operator*(REAL scalar) const
{
  QUAT qn;
  qn.x = x * scalar;
  qn.y = y * scalar;
  qn.z = z * scalar;
  qn.w = w * scalar;
  return qn;
}

/// Multiplies a quaternion by a scalar and stores the result in <b>this</b>
QUAT& QUAT::operator*=(REAL scalar)
{
  x *= scalar;
  y *= scalar;
  z *= scalar;
  w *= scalar;
  return *this;
}

/// Multiplies the 3x3 matrix corresponding to this quaternion by a 3D origin and returns the result in a new 3D origin 
ORIGIN3 QUAT::operator*(const ORIGIN3& o) const
{
  const REAL w2 = w*w;
  const REAL x2 = x*x;
  const REAL y2 = y*y;
  const REAL z2 = z*z;
  const REAL xy = x*y;
  const REAL xz = x*z;
  const REAL yz = y*z;
  const REAL xw = x*w;
  const REAL yw = y*w;
  const REAL zw = z*w;
  return ORIGIN3((-1.0+2.0*(w2+x2))*o.x() + 2.0*((xy-zw)*o.y() + (yw+xz)*o.z()), 2.0*((xy+zw)*o.x() + (-xw+yz)*o.z()) + (-1.0+2.0*(w2+y2))*o.y(), 2.0*((-yw+xz)*o.x() + (xw+yz)*o.y()) + (-1.0+2.0*(w2+z2))*o.z());
}

/// Multiplies <b>this</b> by q and stores the result in <b>this</b>
QUAT& QUAT::operator*=(const QUAT& q)
{
  *this = *this * q;
  return *this;
}

/// Returns the l-infinity norm of the quaternion components
REAL QUAT::norm_inf() const
{
  return std::max(std::fabs(x), std::max(std::fabs(y), std::max(std::fabs(z),
                  std::fabs(w))));
}

/// Calculates the squared magnitude of a quaternion
REAL QUAT::norm_sq() const
{
  return w*w + x*x + y*y + z*z;
}

/// Calculates the magnitude of a quaternion
REAL QUAT::magnitude() const
{
  return safe_sqrt(w*w + x*x + y*y + z*z);
}

/// Computes the product between a quaternion and a scalar
QUAT Ravelin::operator*(REAL scalar, const QUAT& q)
{
  return q * scalar;
}

/// Sends the quaternion to the output stream
std::ostream& Ravelin::operator<<(std::ostream& out, const QUAT& q) 
{
  out << "<" << q.x << ", " << q.y << ", " << q.z << "> " << q.w;

  return out;
}

/// Converts a rotation matrix to a unit Quaternion
QUAT& QUAT::operator=(const MATRIX3& m)
{
  // core computation
  QUAT& q = *this;
  q.w = std::sqrt(std::max((REAL) 0.0, (REAL) 1.0 + m.xx() + m.yy() + m.zz())) * (REAL) 0.5;
  q.x = std::sqrt(std::max((REAL) 0.0, (REAL) 1.0 + m.xx() - m.yy() - m.zz())) * (REAL) 0.5;
  q.y = std::sqrt(std::max((REAL) 0.0, (REAL) 1.0 - m.xx() + m.yy() - m.zz())) * (REAL) 0.5;
  q.z = std::sqrt(std::max((REAL) 0.0, (REAL) 1.0 - m.xx() - m.yy() + m.zz())) * (REAL) 0.5;

  // sign computation
  if (m.zy() - m.yz() < (REAL) 0.0)
    q.x = -q.x;
  if (m.xz() - m.zx() < (REAL) 0.0)
    q.y = -q.y;
  if (m.yx() - m.xy() < (REAL) 0.0)
    q.z = -q.z;

  #ifndef NDEBUG 
  if (!q.unit())
    std::cerr << "QUAT::set() warning!  not a unit quaternion!" << std::endl;
  #endif

  return q;
}

/// Sets quaternion to that represented by an axis-angle representation
QUAT& QUAT::operator=(const AANGLE& a)
{
  const REAL half = a.angle*(REAL) 0.5;
  REAL sina = std::sin(half);
  QUAT& q = *this;
  q.x = a.x * sina;
  q.y = a.y * sina;
  q.z = a.z * sina;
  q.w = std::cos(half);

  #ifndef NDEBUG
  if (!q.unit())
    std::cerr << "QUAT::set() warning!  not a unit quaternion!" << std::endl;
  #endif

  return q;
}


