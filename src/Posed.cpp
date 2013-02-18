/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <Ravelin/cblas.h>
#include <cstring>
#include <cmath>
#include <iomanip>
#include <Ravelin/Constants.h>
#include <Ravelin/Quatd.h>
#include <Ravelin/AAngle.h>
#include <Ravelin/Posed.h>

using namespace Ravelin;

/// Default constructor
/**
 * Sets matrix to the identity matrix
 */
Posed::Posed()
{
  set_identity();
}

/// Constructs a 4x4 homogeneous transformation matrix from a unit quaternion and translation vector
Posed::Posed(const Quatd& q, const Vector3d& v)
{
  Quatd qn = Quatd::normalize(*q);
  set(qn, v);
}

/// Constructs a 4x4 homogeneous transformation matrix from a unit quaternion (for rotation) and zero translation
Posed::Posed(const Quatd& q)
{
  Quatd qn = Quatd::normalize(*q);
  set(qn, ZEROS_3);
}

/// Constructs a 4x4 homogeneous transformation matrix from a rotation matrix and translation vector
Posed::Posed(const Matrix3d& r, const Vector3d& v) 
{
  set(r, v);
}

/// Constructs a 4x4 homogeneous transformation matrix from a rotation matrix and zero translation
Posed::Posed(const Matrix3d& r)
{
  set(r, ZEROS_3d);
}

/// Constructs a 4x4 homogeneous transformation matrix from a axis-angle representation and a translation vector
Posed::Posed(const AAngled& a, const Vector3d& v)
{
  set(a, v);
}

/// Constructs a 4x4 homogeneous transformation matrix from a axis-angle representation (for rotation) and zero translation
Posed::Posed(const AAngled& a)
{
  set(a, ZEROS_3d);
}

/// Constructs a 4x4 homogeneous transformation matrix using identity orientation and a translation vector
Posed::Posed(const Vector3d& v)
{
  set_rotation(Quatd::identity(), v);
}

/// Constructs a matrix from an array
/**
 * \param array an array of 12 Real values in row-major format (all entries
 *        should be given except the last row
 */
Posed::Posed(const Real* array)
{
  set(array);
}

/// Constructs a matrix from 12 values
/**
 * The resulting matrix will appear as follows: <br />
 * m00 m01 m02 m03<br />
 * m10 m11 m12 m13<br />
 * m20 m21 m22 m23<br />
 * 0   0   0   1<br />
 */
Posed::Posed(Real m00, Real m01, Real m02, Real m03, Real m10, Real m11, Real m12, Real m13, Real m20, Real m21, Real m22, Real m23)
{
  const unsigned X = 0, Y = 1, Z = 2, W = 3;

  // set rotation component
  operator()(X,X) = m00;  operator()(X,Y) = m01;  operator()(X,Z) = m02;
  operator()(Y,X) = m10;  operator()(Y,Y) = m11;  operator()(Y,Z) = m12;
  operator()(Z,X) = m20;  operator()(Z,Y) = m21;  operator()(Z,Z) = m22;

  // set translation component
  operator()(X,W) = m03;  operator()(Y,W) = m13;  operator()(Z,W) = m23;
}

/// Interpolates between two 4x4 transforms using spherical linear interpolation
/**
 * \param T1 the matrix to use when t=0
 * \param T2 the matrix to use when t=1
 * \param t a real value in the interval [0,1]
 * \return the interpolated transform
 */
Posed Posed::interpolate(const Posed& T1, const Posed& T2, Real t)
{
  // interpolate the positions
  Vector3d x = T1.x*(1-t) + T2.x*t;

  // interpolate the rotations
  Quatd q = Quatd::slerp(q1.q, q2.q, t);

  // return the new matrix
  return Posed(q, x);
}

/// Sets the matrix using an array of Real values
/**
 * \param array an array of 12 real values in row-major format (fourth row is 
 *        not modified)
 */
void Posed::set(const Real* array)
{
  const unsigned SZ = 3;
  CBLAS::copy(SZ, array, 1, _data, 1);
  CBLAS::copy(SZ, array+SZ, 1, _data+SZ, 1);
  CBLAS::copy(SZ, array+SZ+SZ, 1, _data+SZ+SZ, 1);
  CBLAS::copy(SZ, array+SZ+SZ+SZ, 1, _data+SZ+SZ+SZ, 1);
}

// Sets the matrix to be a 4x4 homogeneous transform from an axis-angle representation and a translation vector
/**
 * \todo make this marginally faster by setting the matrix components directly (eliminates redundancy in setting non-rotation and non-translational components)
 */
void Posed::set(const AAngle* a, const Vector3d*  v)
{
  set_rotation(a);
  set_translation(*v);
}

/// Sets the matrix to be a 4x4 homogeneous transform from a rotation matrix and translation vector
/**
 * \todo make this marginally faster by setting the matrix components directly (eliminates redundancy in setting non-rotation and non-translational components)
 */
void Posed::set(const Matrix3d* m, const Vector3d*  v)
{
  set_rotation(m);
  set_translation(*v);
}

/// Sets the matrix to be a 4x4 homogeneous transform from a unit quaternion and translation vector
/**
 * \todo make this marginally faster by setting the matrix components directly (eliminates redundancy in setting non-rotation and non-translational components)
 */
void Posed::set(const Quatd* q, const Vector3d*  v)
{
  set_rotation(q);
  set_translation(*v);
}

/// Sets the rotation component of a 4x4 homogeneous transform from the given axis-angle representation
void Posed::set_rotation(const AAngle* a)
{
  const unsigned X = 0, Y = 1, Z = 2;
  Real x = a->x;
  Real y = a->y;
  Real z = a->z;
  Real ca = std::cos(a->angle);
  Real sa = std::sin(a->angle);

  // va should be 1 - cos(a); we'll use a slower formula to prevent
  // cancellation error [Ericscon, 2005 (p. 441)]
  const Real SOMEWHAT_NEAR_ZERO = 1e-2;
  Real va = (std::fabs(a->angle) > SOMEWHAT_NEAR_ZERO) ? 1 - ca : (sa*sa)/(1+ca);

  operator()(X,X) = x*x*va + ca;
  operator()(X,Y) = x*y*va - z*sa;
  operator()(X,Z) = x*z*va + y*sa;
  operator()(Y,X) = x*y*va + z*sa;
  operator()(Y,Y) = y*y*va + ca;
  operator()(Y,Z) = y*z*va - x*sa;
  operator()(Z,X) = x*z*va - y*sa;
  operator()(Z,Y) = y*z*va + x*sa;
  operator()(Z,Z) = z*z*va + ca;
}

/// Sets the rotation component of a 4x4 homogeneous transform from the given quaternion 
void Posed::set_rotation(const Quatd* q)
{
  const unsigned X = 0, Y = 1, Z = 2;
  Quatd qn = Quatd::normalize(*q);

  operator()(X,X) = 2*qn.x*qn.x + 2*qn.w*qn.w - 1;
  operator()(X,Y) = 2*qn.x*qn.y - 2*qn.z*qn.w;
  operator()(X,Z) = 2*qn.x*qn.z + 2*qn.y*qn.w;
  operator()(Y,X) = 2*qn.x*qn.y + 2*qn.z*qn.w;
  operator()(Y,Y) = 2*qn.y*qn.y + 2*qn.w*qn.w - 1;
  operator()(Y,Z) = 2*qn.y*qn.z - 2*qn.x*qn.w;
  operator()(Z,X) = 2*qn.x*qn.z - 2*qn.y*qn.w;
  operator()(Z,Y) = 2*qn.y*qn.z + 2*qn.x*qn.w;
  operator()(Z,Z) = 2*qn.z*qn.z + 2*qn.w*qn.w - 1;
}

/// Sets the rotation component of a 4x4 homogeneous transform from the given rotation matrix
void Posed::set_rotation(const Matrix3d* m)
{
  const unsigned SZ = 9;
  for (unsigned i=0; i< SZ; i++)
    _data[i] = m->begin()[i];
}

/// Sets the translation component of a 4x4 homogeneous transform
void Posed::set_translation(Real x, Real y, Real z)
{
  const unsigned TX = 9, TY = 10, TZ = 11;
  _data[TX] = x;
  _data[TY] = y;
  _data[TZ] = z;
}

/// Sets the translation component of a 4x4 homogeneous transform from the given 3-D array
void Posed::set_translation(const Real* array)
{
  const unsigned X = 0, Y = 1, Z = 2, TX = 9, TY = 10, TZ = 11;
  _data[TX] = array[X];
  _data[TY] = array[Y];
  _data[TZ] = array[Z];
}

/// Sets the translation componet of a 4x4 homogeneous transform from the given 3-D vector
void Posed::set_translation(const Vector3d& v)
{
  const unsigned X = 0, Y = 1, Z = 2, TX = 9, TY = 10, TZ = 11;
  _data[TX] = v[X];
  _data[TY] = v[Y];
  _data[TZ] = v[Z];
}

/// Sets this matrix to identity
void Posed::set_identity()
{
  const unsigned XX = 0, YY = 4, ZZ = 8, N = 12;
  for (unsigned i=0; i< N; i++)
    _data[i] = (Real) 0.0;
  _data[XX] = _data[YY] = _data[ZZ] = 1.0;
}

/// Multiplies the 3x3 rotational component of this matrix by a vector (i.e., translation is not added)
Vector3d Posed::mult_vector(const Vector3d& v) const
{
  // note that we artificially set the rows to 3, b/c that is all the data
  // that we have stored
  const unsigned ROWS = 3;

  // this is just regular 3x3 matrix multiplication
  Vector3d result;
  CBLAS::gemv(CblasColMajor, CblasNoTrans, ROWS, ROWS, (Real) 1.0, begin(), ROWS, v.begin(), 1, (Real) 0.0, result.begin(), 1);

  return result;
}

/// Multiplies the transpose of the 3x3 rotational component of this matrix by a vector (i.e., translation is not added)
Vector3d Posed::transpose_mult_vector(const Vector3d& v) const
{
  // note that we artificially set the rows to 3, b/c that is all the data
  // that we have stored
  const unsigned ROWS = 3, X = 0, Y = 1, Z = 2;

  // setup starting rows indices
  const unsigned RX = 0, RY = 3, RZ = 6;

  // this is 3x3 matrix multiplication
  Vector3d result;
  result[X] = CBLAS::dot(ROWS, _data+RX, 1, v.data(), 1);
  result[Y] = CBLAS::dot(ROWS, _data+RY, 1, v.data(), 1);
  result[Z] = CBLAS::dot(ROWS, _data+RZ, 1, v.data(), 1); 

  return result;
}

/// Multiplies a 4x4 matrix by a point 
Vector3d Posed::mult_point(const Vector3d& v) const
{
  // note that we artificially set the rows to 3, b/c that is all the data
  // that we have stored
  const unsigned ROWS = 3;
  const unsigned X = 0, Y = 1, Z = 2, TX = 9, TY = 10, TZ = 11;

  // matrix/vec multiplication is just regular 3x3 matrix multiplication
  Vector3d result;
  CBLAS::gemv(CblasColMajor, CblasNoTrans, ROWS, ROWS, (Real) 1.0, begin(), ROWS, v.begin(), 1, (Real) 0.0, result.begin(), 1);

  // add in the translational components
  result[X] += _data[TX];
  result[Y] += _data[TY];
  result[Z] += _data[TZ];

  return result;
}

/// Multiplies the inverse of the 4x4 matrix by a point
/**
 * \note if multiple inverse_mult_point() operations are performed, it may be
 *       faster to use the inverse_transform() function, followed by 
 *       mult_point().
 */ 
Vector3d Posed::inverse_mult_point(const Vector3d& v) const
{
  // note that we artificially set the rows to 3, b/c that is all the data
  // that we have stored
  const unsigned ROWS = 3, X = 0, Y = 1, Z = 2;

  // setup starting rows indices
  const unsigned RX = 0, RY = 3, RZ = 6, RW = 9;

  // determine R' * v 
  Vector3d result1;
  result1[X] = CBLAS::dot(ROWS, _data+RX, 1, v.data(), 1);
  result1[Y] = CBLAS::dot(ROWS, _data+RY, 1, v.data(), 1);
  result1[Z] = CBLAS::dot(ROWS, _data+RZ, 1, v.data(), 1); 

  // determine R' * -t
  Vector3d result2;
  result2[X] = -CBLAS::dot(ROWS, _data+RX, 1, _data + RW, 1);
  result2[Y] = -CBLAS::dot(ROWS, _data+RY, 1, _data + RW, 1);
  result2[Z] = -CBLAS::dot(ROWS, _data+RZ, 1, _data + RW, 1); 

  // add the two together
  return result1 + result2;
}

/// Special method for inverting a 4x4 transformation matrix in place
void Posed::invert_transform()
{
  // verify that it is a valid transform
  assert(valid_transform(*this));

  // get the new rotation matrix
  q.invert();

  // determine the new translation
  x = q.mult(-x);
}

/// Special method for inverseing a 4x4 transformation matrix
Posed Posed::inverse_transform(const Posed& T)
{
  Posed m(T);
  m.invert_transform();
  return m;
}

/// Multiplies this pose by another
Posed Posed::operator*(const Posed& m) const
{
  // determine the rotation part of the matrix
  const unsigned SZ = 3, XLAT_BEGIN = 9, TX = 9, TY = 10, TZ = 11;
  Posed result;
  CBLAS::gemm(CblasColMajor, CblasNoTrans, CblasNoTrans, SZ, SZ, SZ, (Real) 1.0, _data, SZ, m._data, SZ, (Real) 0.0, result._data, SZ);

  // determine the translation part of the matrix
  CBLAS::gemv(CblasColMajor, CblasNoTrans, SZ, SZ, (Real) 1.0, _data, SZ, m._data + XLAT_BEGIN, 1, (Real) 0.0, result._data + XLAT_BEGIN, 1);
  result._data[TX] += _data[TX];
  result._data[TY] += _data[TY];
  result._data[TZ] += _data[TZ];
  return result;
}

/// Outputs this matrix to the stream
std::ostream& Ravelin::operator<<(std::ostream& out, const Posed& m)
{
  const unsigned ROWS = 3, COLUMNS = 4; 
  for (unsigned i=0; i< ROWS; i++)
  {
    for (unsigned j=0; j< COLUMNS-1; j++)
      out << std::setprecision(15) << m(i,j) << " ";
    out << std::setprecision(15) << m(i,COLUMNS-1) << std::endl;
  }
  out << "0 0 0 1" << std::endl;
   
  return out;
}

