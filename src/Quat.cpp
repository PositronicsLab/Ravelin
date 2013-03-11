/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
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

/// Copy constructor
QUAT::QUAT(const QUAT& q)
{
  *this = q;
}

/// Sets the quaternion from a rotation matrix
QUAT::QUAT(const AANGLE& a)
{
  *this = a;
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
REAL QUAT::operator[](unsigned i) const
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

/// Computes the 3x4 matrix 'L' times a four dimensional array
VECTOR3 QUAT::L_mult(const QUAT& q) const
{
  VECTOR3 v;
  v.x() = -x*q.x + w*q.y + z*q.z - y*q.w;
  v.y() = -y*q.x - z*q.y + w*q.z + x*q.w;
  v.z() = -z*q.x + y*q.y - x*q.z + w*q.w;
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
VECTOR3 QUAT::G_mult(REAL qw, REAL qx, REAL qy, REAL qz) const
{
  VECTOR3 r;
  r.x() = -x*qw + w*qx - z*qy + y*qz;
  r.y() = -y*qw + z*qx + w*qy - x*qz;
  r[2] = -z*qw - y*qx + x*qy + w*qz;
  return r;
}

/// Multiplies the matrix 'G' (transpose) by a vector
QUAT QUAT::G_transpose_mult(const VECTOR3& v) const
{
  QUAT output;
  output.w = -x*v.x() - y*v.y() - z*v.z();
  output.x = w*v.x() + z*v.y() - y*v.z();
  output.y = -z*v.x() + w*v.y() + x*v.z();
  output.z = y*v.x() - x*v.y() + w*v.z();
  return output;
}

/// Multiplies the matrix 'L' (transpose) by a vector
QUAT QUAT::L_transpose_mult(const VECTOR3& v) const
{
  QUAT output;
  output.w = -x*v.x() - y*v.y() - z*v.z();
  output.x = +w*v.x() - z*v.y() + y*v.z();
  output.y = +z*v.x() + w*v.y() - x*v.z();
  output.z = -y*v.x() + x*v.y() + w*v.z();
  return output;
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

  // find the angle between the two
  REAL theta = std::acos(std::fabs(dot));

  REAL sin1at = std::sin((1-alpha)*theta);
  REAL sinat = std::sin(alpha*theta);
  REAL sint_i = (REAL) 1.0/std::sin(theta);
  QUAT qa(q1.x*sin1at, q1.y*sin1at, q1.z*sin1at, q1.w*sin1at);
  QUAT qb(q2.x*sinat, q2.y*sinat, q2.z*sinat, q2.w*sinat);
  if (use_conj)
    qb = -qb;
  QUAT qc = qa + qb;
  qc.x *= sint_i;
  qc.y *= sint_i;
  qc.z *= sint_i;
  qc.w *= sint_i;
  return qc; 
}

/// Computes the inverse orientation represented by this quaternion in place
QUAT& QUAT::inverse()
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

/// Constructs the conjugate of <b>this</b>
QUAT QUAT::operator-() const
{
  return QUAT::conjugate(*this);
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

/// Multiplies the 3x3 matrix corresponding to this quaternion by a 3D vector and returns the result in a new 3D vector
VECTOR3 QUAT::operator*(const VECTOR3& v) const
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
  return VECTOR3((-1.0+2.0*(w2+x2))*v.x() + 2.0*((xy-zw)*v.y() + (yw+xz)*v.z()), 2.0*((xy+zw)*v.x() + (-xw+yz)*v.z()) + (-1.0+2.0*(w2+y2))*v.y(), 2.0*((-yw+xz)*v.x() + (xw+yz)*v.y()) + (-1.0+2.0*(w2+z2))*v.z());
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


