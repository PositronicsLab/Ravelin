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
  set(QUAT::normalize(q), ZEROS_3);
}

/// Constructs a 4x4 homogeneous transformation matrix from a rotation matrix and translation vector
POSE::POSE(const MATRIX3& r, const VECTOR3& v) 
{
  set(r, v);
}

/// Constructs a 4x4 homogeneous transformation matrix from a rotation matrix and zero translation
POSE::POSE(const MATRIX3& r)
{
  set(r, ZEROS_3);
}

/// Constructs a 4x4 homogeneous transformation matrix from a axis-angle representation and a translation vector
POSE::POSE(const AANGLE& a, const VECTOR3& v)
{
  set(a, v);
}

/// Constructs a 4x4 homogeneous transformation matrix from a axis-angle representation (for rotation) and zero translation
POSE::POSE(const AANGLE& a)
{
  set(a, ZEROS_3);
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
  this->x = ZEROS_3;
  return *this;
}

/// Sets this pose from a axis-angle object only (translation will be zero'd) 
POSE& POSE::set(const AANGLE& a)
{
  this->q = a;
  this->x = ZEROS_3;
  return *this;
}

/// Sets this pose from a 3x3 rotation matrix only (translation will be zero'd) 
POSE& POSE::set(const MATRIX3& m)
{
  this->q = m;
  this->x = ZEROS_3;
  return *this;
}

/// Sets this matrix to identity
POSE& POSE::set_identity()
{
  this->q = QUAT::identity();
  this->x = ZEROS_3;
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
  return POSE(q * p.q, q*p.x + x);
}

/// Outputs this matrix to the stream
std::ostream& Ravelin::operator<<(std::ostream& out, const POSE& m)
{
  out << "q: " << m.q << " " << m.x << std::endl;
   
  return out;
}

