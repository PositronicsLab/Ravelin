/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

// ***********************************************************************
// template specializations BEGIN 
// ***********************************************************************

namespace Ravelin {

/// Sets the specified column of the matrix
/**
 * \param i the 0-index of the column
 * \param o a 3D vector 
 */
template <>
MATRIX3& MATRIX3::set_column(unsigned i, const ORIGIN3& o)
{
  const unsigned SZ = 3;
  assert(i < SZ);

  unsigned st = i*SZ;
  _data[st++] = o.x();
  _data[st++] = o.y();
  _data[st] = o.z();
  return *this;
}

/// Gets the specified column of the matrix
/*
 * \param i the 0-index of the column
 */
template <>
ORIGIN3& MATRIX3::get_column(unsigned i, ORIGIN3& result) const
{
  const unsigned SZ = 3;
  assert(i < SZ);

  unsigned st = i*SZ;
  result.x() = _data[st++];
  result.y() = _data[st++];
  result.z() = _data[st];
  return result;
}

/// Sets the specified row of the matrix
/**
 * \param i the 0-index of the row
 * \param o a 3D vector
 * \note the number of columns of this must be three!
 */
template <>
MATRIX3& MATRIX3::set_row(unsigned i, const ORIGIN3& o)
{
  const unsigned SZ = 3;
  assert(i < SZ);

  unsigned st = i;
  _data[st] = o.x();
  st += SZ;
  _data[st] = o.y(); 
  st += SZ;
  _data[st] = o.z(); 
  return *this;
}

/// Gets the specified row of the matrix
/**
 * \param i the 0-index of the row
 */
template <>
ORIGIN3& MATRIX3::get_row(unsigned i, ORIGIN3& result) const
{
  const unsigned SZ = 3;
  assert(i < SZ);

  unsigned st = i;
  result.x() = _data[st];
  st += SZ;
  result.y() = _data[st];
  st += SZ;
  result.z() = _data[st];
  return result;
}
} // end namespace 

// ***********************************************************************
// template specializations END 
// ***********************************************************************

/// Constructs a matrix from an array
/**
 * \param array an array of 9 REAL values in row-major format
 */
MATRIX3::MATRIX3(const REAL* array)
{
  for (unsigned i=0; i< 9; i++)
    _data[i] = array[i];
}

/// Constructs a matrix from a set of values
MATRIX3::MATRIX3(REAL m00, REAL m01, REAL m02, REAL m10, REAL m11, REAL m12, REAL m20, REAL m21, REAL m22)
{
  _data[0] = m00;
  _data[1] = m10;
  _data[2] = m20;
  _data[3] = m01;
  _data[4] = m11;
  _data[5] = m21;
  _data[6] = m02;
  _data[7] = m12;
  _data[8] = m22;
}

/// Multiplies the transpose of this matrix by a vector and returns the result in a new vector
ORIGIN3 MATRIX3::transpose_mult(const ORIGIN3& o) const
{
  ORIGIN3 result;
  result.x() = _data[0]*o.x() + _data[1]*o.y() + _data[2]*o.z();
  result.y() = _data[3]*o.x() + _data[4]*o.y() + _data[5]*o.z();  
  result.z() = _data[6]*o.x() + _data[7]*o.y() + _data[8]*o.z();  
  return result;
}

/// Multiplies the transpose of this matrix by a matrix and returns the result in a new matrix 
MATRIX3 MATRIX3::transpose_mult(const MATRIX3& m) const
{
  MATRIX3 result;

  // do first column
  result.xx() = _data[0]*m.xx() + _data[1]*m.yx() + _data[2]*m.zx();
  result.yx() = _data[3]*m.xx() + _data[4]*m.yx() + _data[5]*m.zx();  
  result.zx() = _data[6]*m.xx() + _data[7]*m.yx() + _data[8]*m.zx();  

  // do second column
  result.xy() = _data[0]*m.xy() + _data[1]*m.yy() + _data[2]*m.zy();
  result.yy() = _data[3]*m.xy() + _data[4]*m.yy() + _data[5]*m.zy();  
  result.zy() = _data[6]*m.xy() + _data[7]*m.yy() + _data[8]*m.zy();  

  // do third column
  result.xz() = _data[0]*m.xz() + _data[1]*m.yz() + _data[2]*m.zz();
  result.yz() = _data[3]*m.xz() + _data[4]*m.yz() + _data[5]*m.zz();  
  result.zz() = _data[6]*m.xz() + _data[7]*m.yz() + _data[8]*m.zz();  

  return result;
}

/// Multiplies this matrix by a vector and returns the result in a new vector
ORIGIN3 MATRIX3::mult(const ORIGIN3& o) const
{
  ORIGIN3 result;
  result.x() = _data[0]*o.x() + _data[3]*o.y() + _data[6]*o.z();
  result.y() = _data[1]*o.x() + _data[4]*o.y() + _data[7]*o.z();  
  result.z() = _data[2]*o.x() + _data[5]*o.y() + _data[8]*o.z();  
  return result;
}

/// Multiplies this matrix by a matrix and returns the result in a new matrix
MATRIX3 MATRIX3::mult(const MATRIX3& m) const
{
  MATRIX3 result;

  // setup first column
  result.xx() = _data[0]*m.xx() + _data[3]*m.yx() + _data[6]*m.zx();
  result.yx() = _data[1]*m.xx() + _data[4]*m.yx() + _data[7]*m.zx();  
  result.zx() = _data[2]*m.xx() + _data[5]*m.yx() + _data[8]*m.zx();  

  // setup second column
  result.xy() = _data[0]*m.xy() + _data[3]*m.yy() + _data[6]*m.zy();
  result.yy() = _data[1]*m.xy() + _data[4]*m.yy() + _data[7]*m.zy();  
  result.zy() = _data[2]*m.xy() + _data[5]*m.yy() + _data[8]*m.zy();  

  // setup third column
  result.xz() = _data[0]*m.xz() + _data[3]*m.yz() + _data[6]*m.zz();
  result.yz() = _data[1]*m.xz() + _data[4]*m.yz() + _data[7]*m.zz();  
  result.zz() = _data[2]*m.xz() + _data[5]*m.yz() + _data[8]*m.zz();  

  return result;
} 

/// Multiplies this matrix by a transposed matrix and returns the result in a new matrix
MATRIX3 MATRIX3::mult_transpose(const MATRIX3& m) const
{
  MATRIX3 result;

  // setup first column
  result.xx() = _data[0]*m.xx() + _data[3]*m.xy() + _data[6]*m.xz();
  result.yx() = _data[1]*m.xx() + _data[4]*m.xy() + _data[7]*m.xz();  
  result.zx() = _data[2]*m.xx() + _data[5]*m.xy() + _data[8]*m.xz();  

  // setup second column
  result.xy() = _data[0]*m.yx() + _data[3]*m.yy() + _data[6]*m.yz();
  result.yy() = _data[1]*m.yx() + _data[4]*m.yy() + _data[7]*m.yz();  
  result.zy() = _data[2]*m.yx() + _data[5]*m.yy() + _data[8]*m.yz();  

  // setup third column
  result.xz() = _data[0]*m.zx() + _data[3]*m.zy() + _data[6]*m.zz();
  result.yz() = _data[1]*m.zx() + _data[4]*m.zy() + _data[7]*m.zz();  
  result.zz() = _data[2]*m.zx() + _data[5]*m.zy() + _data[8]*m.zz();  

  return result;
}

/// Multiplies the transpose of this matrix by a transposed matrix and returns the result in a new matrix
MATRIX3 MATRIX3::transpose_mult_transpose(const MATRIX3& m) const
{
  MATRIX3 result;

  // setup first column
  result.xx() = _data[0]*m.xx() + _data[1]*m.xy() + _data[2]*m.xz();
  result.yx() = _data[3]*m.xx() + _data[4]*m.xy() + _data[5]*m.xz();  
  result.zx() = _data[6]*m.xx() + _data[7]*m.xy() + _data[8]*m.xz();  

  // setup second column
  result.xy() = _data[0]*m.yx() + _data[1]*m.yy() + _data[2]*m.yz();
  result.yy() = _data[3]*m.yx() + _data[4]*m.yy() + _data[5]*m.yz();  
  result.zy() = _data[6]*m.yx() + _data[7]*m.yy() + _data[8]*m.yz();  

  // setup third column
  result.xz() = _data[0]*m.zx() + _data[1]*m.zy() + _data[2]*m.zz();
  result.yz() = _data[3]*m.zx() + _data[4]*m.zy() + _data[5]*m.zz();  
  result.zz() = _data[6]*m.zx() + _data[7]*m.zy() + _data[8]*m.zz();  

  return result;
}

/// Computes the l-infinity norm of this matrix
REAL MATRIX3::norm_inf() const
{
  const unsigned NELMS = 9;

  REAL nrm = (REAL) 0.0;
  for (unsigned i=0; i< NELMS; i++)
    nrm = std::max(nrm, std::fabs(_data[i]));

  return nrm;
}

/// Constructs the vector from a skew-symmetric matrix
VECTOR3 MATRIX3::inverse_skew_symmetric(const MATRIX3& R)
{
  const unsigned X = 0, Y = 1, Z = 2;

  REAL vx = (R.zy() - R.yz()) * (REAL) 0.5;
  REAL vy = (R.xz() - R.zx()) * (REAL) 0.5;
  REAL vz = (R.yx() - R.xy()) * (REAL) 0.5;

  return VECTOR3(vx, vy, vz);
}

/// Constructs a skew-symmetric matrix from the given values
/**
 * The skew symmetric matrix generated will be:
 * | 0   -z   y |
 * | z    0  -x |
 * | -y   x   0 |
 */
MATRIX3 MATRIX3::skew_symmetric(REAL x, REAL y, REAL z)
{
  return MATRIX3(0,-z,y,z,0,-x,-y,x,0);
}

/// Constructs a skew-symmetric matrix from the given values
/**
 * The skew symmetric matrix generated will be:
 * |   0     -v.z()   v.y() |
 * |  v.z()     0    -v.x() |
 * | -v.y()    v.x()    0   |
 */
MATRIX3 MATRIX3::skew_symmetric(const VECTOR3& v)
{
  return skew_symmetric(v.x(), v.y(), v.z());
}

/// Constructs a skew-symmetric matrix from the given values
/**
 * The skew symmetric matrix generated will be:
 * |   0     -v.z()   v.y() |
 * |  v.z()     0    -v.x() |
 * | -v.y()    v.x()    0   |
 */
MATRIX3 MATRIX3::skew_symmetric(const ORIGIN3& v)
{
  return skew_symmetric(v.x(), v.y(), v.z());
}

/// Inverts this matrix 
MATRIX3& MATRIX3::invert()
{
  *this = invert(*this);
  return *this; 
}

/// Calculates the determinant for a 3x3 matrix
REAL MATRIX3::det() const
{
  const unsigned X = 0, Y = 1, Z = 2;
   REAL v1 = xx() * (yy() * zz() - yz() * zy()); 
   REAL v2 = yx() * (xz() * zy() - xy() * zz()); 
   REAL v3 = zx() * (xy() * yz() - xz() * yy()); 
  return v1 + v2 + v3;
}

/// Determines whether this is an orthonormal matrix
bool MATRIX3::is_orthonormal() const
{
  REAL determinant = det();
  return (std::fabs(determinant - (REAL) 1.0) < EPS_FLOAT || std::fabs(determinant + (REAL) 1.0) < EPS_FLOAT);
}

/// Determines the inverse of the given matrix
MATRIX3 MATRIX3::invert(const MATRIX3& m)
{
  // compute the determinant
  REAL determ = 1.0/m.det();
  REAL m00 = determ * (m.yy()*m.zz() - m.zy()*m.yz());
  REAL m01 = determ * (m.zy()*m.xz() - m.xy()*m.zz());
  REAL m02 = determ * (m.xy()*m.yz() - m.yy()*m.xz());
  REAL m10 = determ * (m.yz()*m.zx() - m.yx()*m.zz());
  REAL m11 = determ * (m.xx()*m.zz() - m.zx()*m.xz());
  REAL m12 = determ * (m.yx()*m.xz() - m.xx()*m.yz());
  REAL m20 = determ * (m.yx()*m.zy() - m.zx()*m.yy());
  REAL m21 = determ * (m.zx()*m.xy() - m.xx()*m.zy());
  REAL m22 = determ * (m.xx()*m.yy() - m.xy()*m.yx());
  return MATRIX3(m00, m01, m02, m10, m11, m12, m20, m21, m22);
}

/// Sets this matrix to the rotation matrix of the specified angle around the X axis
void MATRIX3::set_rot_X(REAL angle)
{
  *this = rot_X(angle);
}

/// Returns the rotation matrix of the specified angle around the X axis
MATRIX3 MATRIX3::rot_X(REAL angle)
{
  REAL sina = std::sin(angle);
  REAL cosa = std::cos(angle);
  MATRIX3 r;
  r.xx() = 1;
  r.xy() = 0;
  r.xz() = 0;
  r.yx() = 0;
  r.yy() = cosa;
  r.yz() = -sina;
  r.zx() = 0;
  r.zy() = sina;
  r.zz() = cosa;
  return r;
}

/// Sets this matrix to the rotation matrix of the specified angle around the Y axis
void MATRIX3::set_rot_Y(REAL angle)
{
  *this = rot_Y(angle);
}

/// Returns the rotation matrix of the specified angle around the Y axis
MATRIX3 MATRIX3::rot_Y(REAL angle)
{
  REAL sina = std::sin(angle);
  REAL cosa = std::cos(angle);
  MATRIX3 r;
  r.xx() = cosa;
  r.xy() = 0;
  r.xz() = sina;
  r.yx() = 0;
  r.yy() = 1;
  r.yz() = 0;
  r.zx() = -sina;
  r.zy() = 0;
  r.zz() = cosa;
  return r;
}

/// Sets this matrix to the rotation matrix of the specified angle around the Z axis
void MATRIX3::set_rot_Z(REAL angle)
{
  *this = rot_Z(angle);
}

/// Returns the rotation matrix of the specified angle around the Z axis
MATRIX3 MATRIX3::rot_Z(REAL angle)
{
  REAL sina = std::sin(angle);
  REAL cosa = std::cos(angle);
  MATRIX3 r;
  r.xx() = cosa;
  r.xy() = -sina;
  r.xz() = 0;
  r.yx() = sina;
  r.yy() = cosa;
  r.yz() = 0;
  r.zx() = 0;
  r.zy() = 0;
  r.zz() = 1;
  return r;
}

/// Checks whether this matrix represents a valid rotation and scale for a right-handed coordinate system
bool MATRIX3::valid_rotation_scale(const MATRIX3& R)
{
  const unsigned X = 0, Y = 1, Z = 2;
  const REAL TOLERANCE = std::sqrt(EPS_FLOAT);

  // create vectors for each column
  VECTOR3 cx(R.get_column(X), boost::shared_ptr<const POSE3>());
  VECTOR3 cy(R.get_column(Y), boost::shared_ptr<const POSE3>());
  VECTOR3 cz(R.get_column(Z), boost::shared_ptr<const POSE3>());

  // get the norms for each column
  REAL cx_norm = cx.norm();
  REAL cy_norm = cy.norm();
  REAL cz_norm = cz.norm();

  // check cross product for right-handedness
  VECTOR3 cp = VECTOR3::cross(cx/cx_norm, cy/cy_norm);
  REAL cp_norm = cp.norm();

  return (std::fabs(VECTOR3::dot(cp/cp_norm, cz/cz_norm) - (REAL) 1.0) < TOLERANCE);
}

/// Checks whether this matrix represents a valid rotation for a right-handed coordinate system
bool MATRIX3::valid_rotation(const MATRIX3& R)
{
  REAL TOLERANCE = std::sqrt(EPS_FLOAT);

  // check the determinant of the matrix 
  return (std::fabs(R.det() - 1.0) < TOLERANCE);
}

/// Graham-Schmidt orthonormalization
bool MATRIX3::orthonormalize(VECTOR3& a, VECTOR3& b, VECTOR3& c)
{
  // Gram-Schmidt orthonormalization produces vectors u0, u1, and u2 as
  // follows:
  //
  //   u0 = v0/|v0|
  //   u1 = (v1 - Dot(u0,v1)*u0)/|v1 - Dot(u0,v1)*u0|
  //   u2 = (v2 - Dot(u0,v2)*u0 - Dot(u1,v2)*u1) /
  //          |v2 - Dot(u0,v2)*u0 - Dot(u1,v2)*u1|

  // Compute u0.
  REAL flen = a.norm();
  if (flen == 0.0)
    return false;
  a /= flen;

  // Compute u1.
  REAL fd0 = VECTOR3::dot(a, b);
  b -= fd0*a;
  flen = b.norm();
  if (flen == 0.0)
    return false;
  b /= flen;

  // Compute u2.
  REAL fd1 = VECTOR3::dot(b, c);
  fd0 = VECTOR3::dot(a,c);
  c -= fd0*a + fd1*b;
  flen = c.norm();
  if (flen == 0.0)
    return false;

  c /= flen;
  return true;
}

/// Makes the matrix orthonormal using Gram-Schmidt orthogonalization
bool MATRIX3::orthonormalize()
{
  const unsigned X = 0, Y = 1, Z = 2;

  // The algorithm uses Gram-Schmidt orthogonalization applied to the
  // columns of M.
  VECTOR3 a(xx(),xy(),xz());
  VECTOR3 b(yx(),yy(),yz());
  VECTOR3 c(zx(),zy(),zz());
  if (orthonormalize(a,b,c))
  {
    xx() = a.x();
    xy() = a.y();
    xz() = a.z();
    yx() = b.x();
    yy() = b.y();
    yz() = b.z();
    zx() = c.x();
    zy() = c.y();
    zz() = c.z();
    return true;
  }
  
  return false;
}

/// Gets the specified row of the matrix
/**
 * \param i the 0-index of the row
 *
 */
ORIGIN3 MATRIX3::get_row(unsigned i) const
{
  ORIGIN3 o;
  get_row(i, o);
  return o;
}

/// Gets the specified column of the matrix
/**
 * \param i the 0-index of the column
 */
ORIGIN3 MATRIX3::get_column(unsigned i) const
{
  ORIGIN3 o;
  get_column(i, o);
  return o;
}

/// Sets this matrix to its transpose
void MATRIX3::transpose()
{
  const unsigned YX = 1, ZX = 2, XY = 3, ZY = 5, XZ = 6, YZ = 7;
  std::swap(_data[XY], _data[YX]);
  std::swap(_data[XZ], _data[ZX]);
  std::swap(_data[YZ], _data[ZY]);
}

/// Determines the transpose of this matrix
MATRIX3 MATRIX3::transpose(const MATRIX3& m)
{
  MATRIX3 n = m;
  n.transpose();
  return n;
}

/// Copies a matrix to this one
MATRIX3& MATRIX3::operator=(const MATRIX3& m)
{
  const unsigned SZ = 9;
  for (unsigned i=0; i< SZ; i++)
    _data[i] = m._data[i];
  return *this;
}

/// Copies a matrix to this one
MATRIX3& MATRIX3::operator=(const MATRIXN& m)
{
  #ifndef NEXCEPT
  if (m.rows() != 3 || m.columns() != 3)
    throw MissizeException();
  #endif
  const REAL* data = m.data();
  const unsigned SZ = 9;
  for (unsigned i=0; i< SZ; i++)
    _data[i] = data[i];
  return *this;
}

/// Copies a matrix to this one
MATRIX3& MATRIX3::operator=(const SHAREDMATRIXN& m)
{
  #ifndef NEXCEPT
  if (m.rows() != 3 || m.columns() != 3)
    throw MissizeException();
  #endif
  const unsigned SZ = 9;
  CONST_COLUMN_ITERATOR miter = m.column_iterator_begin();
  for (unsigned i=0; i< SZ; i++)
    _data[i] = *miter++;
  return *this;
}

/// Copies a matrix to this one
MATRIX3& MATRIX3::operator=(const CONST_SHAREDMATRIXN& m)
{
  #ifndef NEXCEPT
  if (m.rows() != 3 || m.columns() != 3)
    throw MissizeException();
  #endif
  const unsigned SZ = 9;
  CONST_COLUMN_ITERATOR miter = m.column_iterator_begin();
  for (unsigned i=0; i< SZ; i++)
    _data[i] = *miter++;
  return *this;
}

/// Multiplies this matrix by a scalar in place
MATRIX3& MATRIX3::operator*=(REAL scalar)
{
  const unsigned SZ = 3*3;
  CBLAS::scal(SZ, scalar, data(), 1);
  return *this;
}

/// Returns the negation of this matrix
MATRIX3 MATRIX3::operator-() const
{
  MATRIX3 m;
  std::transform(_data, _data+9, m._data, std::negate<REAL>());
  return m;
}

/// Adds m to this in place
MATRIX3& MATRIX3::operator+=(const MATRIX3& m)
{
  std::transform(_data, _data+9, m._data, _data, std::plus<REAL>());
  return *this;
}

/// Subtracts m from this in place
MATRIX3& MATRIX3::operator-=(const MATRIX3& m)
{
  std::transform(_data, _data+9, m._data, _data, std::minus<REAL>());
  return *this;
}

/// Sets this matrix to identity
MATRIX3& MATRIX3::set_identity()
{
  const unsigned XX = 0, YY = 4, ZZ = 8;
  xy() = xz() = yz() = yx() = zx() = zy() = (REAL) 0.0;
  _data[XX] = _data[YY] = _data[ZZ] = (REAL) 1.0;
  return *this;
}

/// Checks whether this matrix is symmetric
bool MATRIX3::is_symmetric(REAL tolerance) const
{
  // make sure tolerance is non-negative
  tolerance = std::fabs(tolerance);

  // check symmetric
  for (unsigned i=0; i< rows()-1; i++)
    for (unsigned j=i+1; j< rows(); j++)
      if (!OPS::rel_equal(operator()(i,j), operator()(j,i), tolerance))
        return false;

  return true;
}

/// Calculates the differential between two rotations
VECTOR3 MATRIX3::calc_differential(const MATRIX3& R1, const MATRIX3& R2)
{
  const unsigned X = 0, Y = 1, Z = 2;
  VECTOR3 R1x(R1.get_column(X), boost::shared_ptr<POSE3>());
  VECTOR3 R1y(R1.get_column(Y), boost::shared_ptr<POSE3>());
  VECTOR3 R1z(R1.get_column(Z), boost::shared_ptr<POSE3>());
  VECTOR3 R2x(R2.get_column(X), boost::shared_ptr<POSE3>());
  VECTOR3 R2y(R2.get_column(Y), boost::shared_ptr<POSE3>());
  VECTOR3 R2z(R2.get_column(Z), boost::shared_ptr<POSE3>());
  return (VECTOR3::cross(R1x,R2x) + VECTOR3::cross(R1y,R2y) + 
          VECTOR3::cross(R1z,R2z))*0.5;
}

/// Outputs this matrix to the stream
std::ostream& Ravelin::operator<<(std::ostream& out, const MATRIX3& m)
{
  const unsigned SZ = 3; 
  for (unsigned i=0; i< SZ; i++)
  {
    for (unsigned j=0; j< SZ-1; j++)
      out << std::setprecision(15) << m(i,j) << " ";
    out << std::setprecision(15) << m(i,SZ-1) << std::endl;
  }
   
  return out;
}

MATRIX3& MATRIX3::operator=(const QUAT& q)
{
  const unsigned X = 0, Y = 1, Z = 2;
 
  // verify that the quaternion is normalized
  assert(std::fabs(q.magnitude()) - 1 < EPS_FLOAT);

  // setup repeated products
  const REAL xx = q.x*q.x;
  const REAL xy = q.x*q.y;
  const REAL xz = q.x*q.z;
  const REAL xw = q.x*q.w;
  const REAL yy = q.y*q.y;
  const REAL yz = q.y*q.z;
  const REAL yw = q.y*q.w;
  const REAL zz = q.z*q.z;
  const REAL zw = q.z*q.w; 
  const REAL ww = q.w*q.w;

  MATRIX3& m = *this;
  m.xx() = 2*(xx + ww) - 1;
  m.xy() = 2*(xy - zw);
  m.xz() = 2*(xz + yw);
  m.yx() = 2*(xy + zw);
  m.yy() = 2*(yy + ww) - 1;
  m.yz() = 2*(yz - xw);
  m.zx() = 2*(xz - yw);
  m.zy() = 2*(yz + xw);
  m.zz() = 2*(zz + ww) - 1;
  return m;
}

MATRIX3& MATRIX3::operator=(const AANGLE& a)
{
  const unsigned X = 0, Y = 1, Z = 2;
  REAL x = a.x;
  REAL y = a.y;
  REAL z = a.z;
  REAL ca = cos(a.angle);
  REAL sa = sin(a.angle);

  // va should be 1 - cos(a); we'll use a slower formula to prevent
  // cancellation error [Ericscon, 2005 (p. 441)]
  const REAL SOMEWHAT_EPS_FLOAT = 1e-2;
  REAL va = (std::fabs(a.angle) > SOMEWHAT_EPS_FLOAT) ? 1 - ca : (sa*sa)/(1+ca);

  // setup the matrix
  MATRIX3& m = *this;
  m.xx() = x*x*va + ca;
  m.xy() = x*y*va - z*sa;
  m.xz() = x*z*va + y*sa;
  m.yx() = x*y*va + z*sa;
  m.yy() = y*y*va + ca;
  m.yz() = y*z*va - x*sa;
  m.zx() = x*z*va - y*sa;
  m.zy() = y*z*va + x*sa;
  m.zz() = z*z*va + ca;
  return m;
}

REAL& MATRIX3::operator()(unsigned i, unsigned j)
{
  #ifndef NEXCEPT
  if (i >= 3 || j >= 3)
    throw InvalidIndexException();
  #endif
  return _data[j*3+i];
}

const REAL& MATRIX3::operator()(unsigned i, unsigned j) const
{
  #ifndef NEXCEPT
  if (i >= 3 || j >= 3)
    throw InvalidIndexException();
  #endif
  return _data[j*3+i];
}

REAL* MATRIX3::data(unsigned i)
{
  #ifndef NEXCEPT
  if (i >= 9)
    throw InvalidIndexException();
  #endif
  return &_data[i];
}

const REAL* MATRIX3::data(unsigned i) const
{
  #ifndef NEXCEPT
  if (i >= 9)
    throw InvalidIndexException();
  #endif
  return &_data[i];
}

MATRIX3& MATRIX3::resize(unsigned m, unsigned n, bool preserve)
{
  #ifndef NEXCEPT
  if (m != 3 || n != 3)
    throw std::runtime_error("Attempt to resize fixed-length vector!");
  #endif

  return *this;
}
 

